import os
import sys

import pandas as pd
import openpyxl
import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, Eq, solve
from enum import Enum

class Type(Enum):
    NM808 = 1
    NM248 = 2


def read_excel_sheet(location) -> (pd.ExcelFile, dict):
    xl = pd.ExcelFile(location)

    sheet_names = xl.sheet_names

    sheets = {}

    for sheet in sheet_names:
        sheets[sheet] = xl.parse(sheet, header=None)

    return xl, sheets


def write_data_to_excel(sheets, output_location) -> None:
    # Create a new Excel writer object
    writer = pd.ExcelWriter(output_location, engine='xlsxwriter')

    # Write each DataFrame to a different sheet in the Excel file
    for sheet_name, data in sheets.items():
        # Increase the precision of the numbers
        data.to_excel(writer, sheet_name=sheet_name)

    # Save the Excel file
    writer.close()


def reset_indexes(*data):
    result = []
    for df in data:
        if isinstance(df, pd.Series):
            df = df.to_frame()
        df = df.reset_index(drop=True)
        if df.shape[1] == 1:
            df.columns = [0]
        else:
            df.columns = range(df.shape[1])
        result.append(df)
    if len(result) == 1:
        return result[0]
    return tuple(result)

def calculate_std_808(row):
    # Apply the operation to specific cells and calculate the standard deviation
    return np.std([(cell * 1000000) / 7900 for cell in row], ddof=0)


def calculate_std_248(row):
    # Apply the operation to specific cells and calculate the standard deviation
    return np.std([((cell + 0.1058) / 4e-5) for cell in row], ddof=0)


def get_slope_intercept(equation):
    equation = equation.replace("^", "**")
    equation = equation.replace("y=", "")

    x, y = symbols('x y')  # Define 'x' and 'y'

    # Parse the equation
    eq = Eq(y, eval(equation))

    # Solve the equation for y
    solution = solve(eq, y)

    # The slope is the coefficient of x
    slope = solution[0].coeff(x)

    # The intercept is the constant term
    intercept = solution[0].subs(x, 0)

    intercept = abs(intercept)

    return float(slope), float(intercept)


def convert_to_numeric_dataframes(exit_if_null, *data) -> tuple[pd.DataFrame, ...]:
    numeric = []
    for each in data:
        each = each.apply(pd.to_numeric, errors='coerce')
        if exit_if_null:
            if each.isnull().values.any():
                print("Please check the data provided. There are blank cells in the data.")
                print("The blank cells are located at the following row and column:")
                print(each[each.isnull().any(axis=1)].index)
                print(each[each.isnull().any(axis=1)].columns)
                sys.exit(1)
        each = each.dropna()
        numeric.append(each)
    return tuple(numeric)

def concentration_calculations(data, loaded, control) -> tuple[float, float]:
    interior_volume = data.iloc[3, 1]
    exterior_volume = data.iloc[4, 1]
    loaded_swnts = (loaded * (interior_volume / exterior_volume))
    control_swnts = (control * (interior_volume / exterior_volume))
    return loaded_swnts, control_swnts

def find_row(data, nanometer) -> slice:
    start_row, end_row = 0, 0
    time_column = data.iloc[:, 3]
    if nanometer == Type.NM248:
        for index, value in enumerate(time_column):
            if value == 'Absorption 248 nm':
                start_row = index + 2
            if value == 'Absorption 808 nm':
                end_row = index - 4
                break
        return slice(start_row, end_row)
    elif nanometer == Type.NM808:
        for index, value in enumerate(time_column):
            if value == 'Absorption 808 nm':
                start_row = index + 2
        end_row = time_column.size
        return slice(start_row, end_row)

def find_columns(data) -> dict[str, list]:
    header_row = data.iloc[1, :17]
    cols = {"Time": [], "R848": [], "SWNTs Control": [], "Drug Control": [], "PBS Control": []}
    for index, value in enumerate(header_row):
        value = str(value)
        if value.startswith("Time"):
            cols["Time"].append(index)
        elif value.startswith("R848"):
            cols["R848"].append(index)
        elif value.startswith("SWNT"):
            cols["SWNTs Control"].append(index)
        elif value.startswith("Drug Con"):
            cols["Drug Control"].append(index)
        elif value.startswith("PBS"):
            cols["PBS Control"].append(index)
    return cols


def calculate(data, replication, condensed) -> pd.DataFrame:
    dfs = []
    factor_loaded_swnts, factor_control = concentration_calculations(data, data.iloc[9, 1], data.iloc[11, 1])
    slope, intercept = get_slope_intercept(data.iloc[6, 1])
    pbs_control_248, pbs_control_808, pbs_calculations = None, None, False

    rows_248 = find_row(data, Type.NM248)
    if rows_248.start == 0 and rows_248.stop == 0:
        print("Invalid data provided. Invalid values found for 248 nm. Check the Excel sheet and try again.")
        sys.exit(1)

    rows_808 = find_row(data, Type.NM808)
    if rows_808.start == 0 and rows_808.stop == 0:
        print("Invalid data provided. Invalid values found for 808 nm. Check the Excel sheet and try again.")
        sys.exit(1)

    cols = find_columns(data)

    if len(cols['PBS Control']) != 0:
        pbs_calculations = True
    time = data.iloc[rows_248.start:rows_248.stop, cols['Time']]
    time, absorption_808_r848, absorption_808_swnt_control, absorption_808_drug_control = (
        reset_indexes(time, data.iloc[rows_808.start:rows_808.stop, cols['R848']],
                      data.iloc[rows_808.start:rows_808.stop, cols['SWNTs Control']],
                      data.iloc[rows_808.start:rows_808.stop, cols['Drug Control']]))
    absorption_248_r848, absorption_248_swnt_control, absorption_248_drug_control = (
        reset_indexes(data.iloc[rows_248.start:rows_248.stop, cols['R848']],
                      data.iloc[rows_248.start:rows_248.stop, cols['SWNTs Control']],
                      data.iloc[rows_248.start:rows_248.stop, cols['Drug Control']]))

    if pbs_calculations:
        pbs_control_248, pbs_control_808 = reset_indexes(data.iloc[rows_248.start:rows_248.stop, cols['PBS Control']],
                                                         data.iloc[rows_808.start:rows_808.stop, cols['PBS Control']])

    percentage_calculations = pd.DataFrame()

    graph_nm = pd.DataFrame()
    graph_perc = pd.DataFrame()

    (absorption_808_r848, absorption_808_swnt_control, absorption_808_drug_control,
     absorption_248_r848, absorption_248_swnt_control, absorption_248_drug_control,
     time) = convert_to_numeric_dataframes(True, absorption_808_r848, absorption_808_swnt_control,
                                           absorption_808_drug_control,
                                           absorption_248_r848, absorption_248_swnt_control,
                                           absorption_248_drug_control, time)

    if pbs_calculations:
        pbs_control_248, pbs_control_808 = convert_to_numeric_dataframes(True, pbs_control_248, pbs_control_808)

    df = pd.DataFrame()
    df['Time (h)'] = time
    graph_nm['Time'] = time
    graph_perc['Time'] = time
    dfs.append(df)

    df = pd.DataFrame()
    df['R848 SWNTs - 808nm'] = ((absorption_808_r848.sum(axis=1) / replication) * 1000000) / 7900
    df['STD (%)'] = absorption_808_r848.apply(calculate_std_808, axis=1)
    if not condensed:
        dfs.append(df)

    df = pd.DataFrame()
    df['SWNTs Control - 808nm'] = ((absorption_808_swnt_control.sum(axis=1) / replication) * 1000000) / 7900
    df['STD (%)'] = absorption_808_swnt_control.apply(calculate_std_808, axis=1)
    if not condensed:
        dfs.append(df)

    df = pd.DataFrame()
    df['Drug Control - 808nm'] = ((absorption_808_drug_control.sum(axis=1) / replication) * 1000000) / 7900
    df['STD(%)'] = absorption_808_drug_control.apply(calculate_std_808, axis=1)
    if not condensed:
        dfs.append(df)

    df = pd.DataFrame()
    df['R848 SWNTs - 248nm'] = ((absorption_248_r848.sum(axis=1) / replication) + intercept) / slope
    percentage_calculations['R848 SWNTS'] = df['R848 SWNTs - 248nm']
    graph_nm['R848 SWNTS'] = df['R848 SWNTs - 248nm']
    df['STD (%)'] = absorption_248_r848.apply(calculate_std_248, axis=1)
    percentage_calculations['R848 STD'] = df['STD (%)']
    dfs.append(df)

    df = pd.DataFrame()
    df['SWNTs Control - 248nm'] = ((absorption_248_swnt_control.sum(axis=1) / replication) + intercept) / slope
    percentage_calculations['SWNTS Control'] = df['SWNTs Control - 248nm']
    graph_nm['SWNTS Control'] = df['SWNTs Control - 248nm']
    df['STD (%)'] = absorption_248_swnt_control.apply(calculate_std_248, axis=1)
    dfs.append(df)

    df = pd.DataFrame()
    df['Drug Control - 248nm'] = ((absorption_248_drug_control.sum(axis=1) / replication) + intercept) / slope
    percentage_calculations['Drug Control'] = df['Drug Control - 248nm']
    graph_nm['Drug Control'] = df['Drug Control - 248nm']
    df['STD (%)'] = absorption_248_drug_control.apply(calculate_std_248, axis=1)
    percentage_calculations['Drug Control STD'] = df['STD (%)']
    dfs.append(df)

    df = pd.DataFrame()
    df['R848 SWNTS - 248nm (%)'] = (((percentage_calculations['R848 SWNTS'] -
                               percentage_calculations["SWNTS Control"]) / factor_loaded_swnts) * 100)
    graph_perc['R848 SWNTs'] = df['R848 SWNTS - 248nm (%)']
    df['STD (%)'] = (percentage_calculations['R848 STD'] / factor_loaded_swnts) * 100
    dfs.append(df)

    df = pd.DataFrame()
    df['Drug Control - 248nm (%)'] = ((percentage_calculations['Drug Control'] - 2645) / factor_control) * 100
    graph_perc['Drug Control'] = df['Drug Control - 248nm (%)']
    df['STD (%)'] = (percentage_calculations['Drug Control STD'] / factor_control) * 100
    dfs.append(df)

    if pbs_calculations:
        df = pd.DataFrame()
        df['PBS Control - 248nm'] = ((pbs_control_248.sum(axis=1) / replication) + intercept) / slope
        percentage_calculations['PBS Control 248'] = df['PBS Control - 248nm']
        graph_nm['PBS Control'] = df['PBS Control - 248nm']
        df['STD'] = pbs_control_248.apply(calculate_std_248, axis=1)
        dfs.append(df)

        df = pd.DataFrame()
        df['PBS Control - 808nm'] = ((pbs_control_808.sum(axis=1) / replication) * 1000000) / 7900
        df['STD (%)'] = pbs_control_808.apply(calculate_std_808, axis=1)
        if not condensed:
            dfs.append(df)

    values = pd.concat(dfs, axis=1)

    graph_nm, graph_perc = convert_to_numeric_dataframes(False,
                                                         graph_nm, graph_perc)

    plt.figure()
    for column in graph_nm.columns:
        if not pbs_calculations and column == 'PBS Control':
            continue
        plt.scatter(graph_nm['Time'], graph_nm[column], label=column)
        z = np.polyfit(graph_nm['Time'], graph_nm[column], 1)
        p = np.poly1d(z)
        plt.plot(graph_nm['Time'], p(graph_nm['Time']), linestyle='--')

    plt.title('R848 Release (nM)')
    plt.xlabel('Time (h)')
    plt.ylabel('Average (nM)')
    plt.legend()
    plt.grid(True)

    plt.savefig('R848 Release (nM).png', bbox_inches='tight')
    plt.close()

    plt.figure()

    for column in graph_perc.columns:
        plt.scatter(graph_perc['Time'], graph_perc[column], label=column)
        z = np.polyfit(graph_perc['Time'], graph_perc[column], 1)
        p = np.poly1d(z)
        plt.plot(graph_perc['Time'], p(graph_perc['Time']), linestyle='--')

    plt.title("R848 Release (%)")
    plt.xlabel('Time (h)')
    plt.ylabel('Average (%)')
    plt.legend()
    plt.grid(True)

    plt.savefig('R848 Release (%).png', bbox_inches='tight')
    plt.close()

    return values


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python resiquimod_calculations.py <excel sheet location(str)> <sheet_name(str)> <condensed(bool)>")
        sys.exit(1)
    excel_sheet_location = sys.argv[1]
    required_sheet_name = sys.argv[2]
    condensed = True
    if sys.argv == 5:
        if sys.argv[4] not in ['True', 'False']:
            print("The fourth argument must be either True or False")
            sys.exit(1)
        condensed = bool(sys.argv[4])
    if not (os.path.exists(excel_sheet_location) and os.path.isfile(excel_sheet_location)):
        if not excel_sheet_location.endswith('.xlsx'):
            print("The file provided is not an Excel file")
            sys.exit(1)
        print("The file provided does not exist")
        sys.exit(1)
    file, full_data = read_excel_sheet(excel_sheet_location)
    if required_sheet_name not in full_data:
        print("The sheet name provided does not exist in the Excel file")
        sys.exit(1)
    data = full_data[required_sheet_name]
    replication = data.iloc[7, 1]

    calculated_output = calculate(data, replication, condensed)

    pH_value = data.iloc[5, 1]

    output_file = f"R848 Release pH={pH_value}.xlsx"

    # Now write the DataFrames to an Excel file
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        # Write each DataFrame to a different sheet or the same sheet in Excel
        calculated_output.to_excel(writer, sheet_name=f"R848 Release pH={pH_value}", startrow=0, startcol=0, index=False)

        # Access the xlsxwriter workbook and worksheet objects
        workbook = writer.book
        worksheet = writer.sheets[f"R848 Release pH={pH_value}"]

        # Calculate the position where you want to insert your plots
        # For example, below the DataFrames
        image_row = len(data) + 3

        # Insert the plots
        worksheet.insert_image(image_row, 0, 'R848 Release (nM).png')
        worksheet.insert_image(image_row, 4, 'R848 Release (%).png')

        for column in range(calculated_output.shape[1]):
            worksheet.set_column(column, column, 24)  # Set the column width to 15

    os.remove('R848 Release (nM).png')
    os.remove('R848 Release (%).png')
