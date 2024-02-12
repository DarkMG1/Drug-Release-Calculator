import os
import sys

import pandas as pd
import openpyxl, xlsxwriter
import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, Eq, solve
from enum import Enum
import importlib.util

class Type(Enum):
    NM808 = 1
    NM248 = 2

def variable_check(var, desired_type):
    if not isinstance(var, desired_type):
        print(f"Invalid value found for {var}. Please check the Excel sheet and try again.")
        sys.exit(1)



def read_excel_sheet(location) -> (pd.ExcelFile, dict):
    xl = pd.ExcelFile(location)

    sheet_names = xl.sheet_names

    sheets = {}

    for sheet in sheet_names:
        sheets[sheet] = xl.parse(sheet, header=None)

    return xl, sheets

def is_module_installed(module_name) -> bool:
    spec = importlib.util.find_spec(module_name)
    return spec is not None

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
        new_data = each.apply(pd.to_numeric, errors='coerce')
        if exit_if_null:
            if new_data.isnull().values.any():
                print("Please check the data provided. There are blank cells in the data.")
                print(new_data)
                print("Above is the data that is causing the issue. Please check the Excel sheet and try again.")
                sys.exit(1)
        new_data = new_data.dropna()
        numeric.append(new_data)
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
                end_row = index - 1
                break
        return slice(start_row, end_row)
    elif nanometer == Type.NM808:
        for index, value in enumerate(time_column):
            if value == 'Absorption 808 nm':
                start_row = index + 2
        end_row = time_column.size
        return slice(start_row, end_row)

def find_columns(data, drug_name) -> dict[str, list]:
    header_row = data.iloc[1, :17]
    cols = {"Time": [], drug_name: [], "SWNTs Control": [], "Drug Control": [], "PBS Control": []}
    for index, value in enumerate(header_row):
        value = str(value)
        if value.startswith("Time"):
            cols["Time"].append(index)
        elif value.startswith(drug_name):
            cols[drug_name].append(index)
        elif value.startswith("SWNT"):
            cols["SWNTs Control"].append(index)
        elif value.startswith("Drug Con"):
            cols["Drug Control"].append(index)
        elif value.startswith("PBS"):
            cols["PBS Control"].append(index)
    return cols


def calculate(data, replication, drug_name) -> tuple[pd.DataFrame, pd.DataFrame]:
    nm_248 = []
    nm_808 = []
    if data.iloc[9, 1] == 0 or data.iloc[11, 1] == 0 or variable_check(data.iloc[9, 1], int) or variable_check(data.iloc[11, 1], int):
        print("Invalid values found for loaded SWNTs and control SWNTs. Please check the Excel sheet and try again.")
        sys.exit(1)
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

    cols = find_columns(data, drug_name)

    if len(cols['PBS Control']) != 0:
        pbs_calculations = True
    time = data.iloc[rows_248.start:rows_248.stop, cols['Time']]
    time, absorption_808_loaded_drug, absorption_808_swnt_control, absorption_808_drug_control = (
        reset_indexes(time, data.iloc[rows_808.start:rows_808.stop, cols[drug_name]],
                      data.iloc[rows_808.start:rows_808.stop, cols['SWNTs Control']],
                      data.iloc[rows_808.start:rows_808.stop, cols['Drug Control']]))
    absorption_248_loaded_drug, absorption_248_swnt_control, absorption_248_drug_control = (
        reset_indexes(data.iloc[rows_248.start:rows_248.stop, cols[drug_name]],
                      data.iloc[rows_248.start:rows_248.stop, cols['SWNTs Control']],
                      data.iloc[rows_248.start:rows_248.stop, cols['Drug Control']]))

    if pbs_calculations:
        pbs_control_248, pbs_control_808 = reset_indexes(data.iloc[rows_248.start:rows_248.stop, cols['PBS Control']],
                                                         data.iloc[rows_808.start:rows_808.stop, cols['PBS Control']])

    percentage_calculations_248 = pd.DataFrame()
    percentage_calculations_808 = pd.DataFrame()

    graph_nm_248 = pd.DataFrame()
    graph_perc_248 = pd.DataFrame()

    graph_nm_808 = pd.DataFrame()
    graph_perc_808 = pd.DataFrame()

    (absorption_808_loaded_drug, absorption_808_swnt_control, absorption_808_drug_control,
     absorption_248_loaded_drug, absorption_248_swnt_control, absorption_248_drug_control,
     time) = convert_to_numeric_dataframes(True, absorption_808_loaded_drug, absorption_808_swnt_control,
                                           absorption_808_drug_control,
                                           absorption_248_loaded_drug, absorption_248_swnt_control,
                                           absorption_248_drug_control, time)

    if pbs_calculations:
        pbs_control_248, pbs_control_808 = convert_to_numeric_dataframes(True, pbs_control_248, pbs_control_808)

    df = pd.DataFrame()
    df['Time (h)'] = time
    graph_nm_248['Time'] = time
    graph_perc_248['Time'] = time
    graph_nm_808['Time'] = time
    graph_perc_808['Time'] = time
    nm_248.append(df)
    nm_808.append(df)

    df = pd.DataFrame()
    df[f'{drug_name} SWNTs - 808nm'] = ((absorption_808_loaded_drug.sum(axis=1) / replication) * 1000000) / 7900
    percentage_calculations_808[drug_name + ' SWNTs'] = df[f'{drug_name} SWNTs - 808nm']
    graph_nm_808[f'{drug_name} SWNTs'] = df[f'{drug_name} SWNTs - 808nm']
    df['STD (%)'] = absorption_808_loaded_drug.apply(calculate_std_808, axis=1)
    percentage_calculations_808[f'{drug_name} STD'] = df['STD (%)']
    nm_808.append(df)

    df = pd.DataFrame()
    df['SWNTs Control - 808nm'] = ((absorption_808_swnt_control.sum(axis=1) / replication) * 1000000) / 7900
    percentage_calculations_808['SWNTs Control'] = df['SWNTs Control - 808nm']
    graph_nm_808['SWNTs Control'] = df['SWNTs Control - 808nm']
    df['STD (%)'] = absorption_808_swnt_control.apply(calculate_std_808, axis=1)
    percentage_calculations_808['SWNTs Control STD'] = df['STD (%)']
    nm_808.append(df)

    df = pd.DataFrame()
    df['Drug Control - 808nm'] = ((absorption_808_drug_control.sum(axis=1) / replication) * 1000000) / 7900
    percentage_calculations_808['Drug Control'] = df['Drug Control - 808nm']
    graph_nm_808['Drug Control'] = df['Drug Control - 808nm']
    df['STD(%)'] = absorption_808_drug_control.apply(calculate_std_808, axis=1)
    percentage_calculations_808['Drug Control STD'] = df['STD(%)']
    nm_808.append(df)

    df = pd.DataFrame()
    df[f'{drug_name} SWNTs - 248nm'] = ((absorption_248_loaded_drug.sum(axis=1) / replication) + intercept) / slope
    percentage_calculations_248[f'{drug_name} SWNTs'] = df[f'{drug_name} SWNTs - 248nm']
    graph_nm_248[f'{drug_name} SWNTs'] = df[f'{drug_name} SWNTs - 248nm']
    df['STD (%)'] = absorption_248_loaded_drug.apply(calculate_std_248, axis=1)
    percentage_calculations_248[f'{drug_name} STD'] = df['STD (%)']
    nm_248.append(df)

    df = pd.DataFrame()
    df['SWNTs Control - 248nm'] = ((absorption_248_swnt_control.sum(axis=1) / replication) + intercept) / slope
    percentage_calculations_248['SWNTs Control'] = df['SWNTs Control - 248nm']
    graph_nm_248['SWNTs Control'] = df['SWNTs Control - 248nm']
    df['STD (%)'] = absorption_248_swnt_control.apply(calculate_std_248, axis=1)
    nm_248.append(df)

    df = pd.DataFrame()
    df['Drug Control - 248nm'] = ((absorption_248_drug_control.sum(axis=1) / replication) + intercept) / slope
    percentage_calculations_248['Drug Control'] = df['Drug Control - 248nm']
    graph_nm_248['Drug Control'] = df['Drug Control - 248nm']
    df['STD (%)'] = absorption_248_drug_control.apply(calculate_std_248, axis=1)
    percentage_calculations_248['Drug Control STD'] = df['STD (%)']
    nm_248.append(df)

    if pbs_calculations:
        df = pd.DataFrame()
        df['PBS Control - 248nm'] = ((pbs_control_248.sum(axis=1) / replication) + intercept) / slope
        percentage_calculations_248['PBS Control 248'] = df['PBS Control - 248nm']
        graph_nm_248['PBS Control'] = df['PBS Control - 248nm']
        df['STD'] = pbs_control_248.apply(calculate_std_248, axis=1)
        nm_248.append(df)

        df = pd.DataFrame()
        df['PBS Control - 808nm'] = ((pbs_control_808.sum(axis=1) / replication) * 1000000) / 7900
        percentage_calculations_808['PBS Control 808'] = df['PBS Control - 808nm']
        df['STD (%)'] = pbs_control_808.apply(calculate_std_808, axis=1)
        nm_808.append(df)

    df = pd.DataFrame()
    df[f'{drug_name} SWNTs - 248nm (%)'] = (((percentage_calculations_248[f'{drug_name} SWNTs'] -
                               percentage_calculations_248["SWNTs Control"]) / factor_loaded_swnts) * 100)
    graph_perc_248[f'{drug_name} SWNTs'] = df[f'{drug_name} SWNTs - 248nm (%)']
    df['STD (%)'] = (percentage_calculations_248[f'{drug_name} STD'] / factor_loaded_swnts) * 100
    nm_248.append(df)

    df = pd.DataFrame()
    df[f'{drug_name} SWNTs - 808nm (%)'] = (((percentage_calculations_808[f'{drug_name} SWNTs'] -
                               percentage_calculations_808["SWNTs Control"]) / factor_loaded_swnts) * 100)
    graph_perc_808[f'{drug_name} SWNTs'] = df[f'{drug_name} SWNTs - 808nm (%)']
    df['STD (%)'] = (percentage_calculations_808[f'{drug_name} STD'] / factor_loaded_swnts) * 100
    nm_808.append(df)



    if pbs_calculations:
        df = pd.DataFrame()
        df['Drug Control - 248nm (%)'] = ((percentage_calculations_248['Drug Control'] - percentage_calculations_248['PBS Control 248'])
                                          / factor_control) * 100
        graph_perc_248['Drug Control'] = df['Drug Control - 248nm (%)']
        df['STD (%)'] = (percentage_calculations_248['Drug Control STD'] / factor_control) * 100
        nm_248.append(df)

        df = pd.DataFrame()
        df['Drug Control - 808nm (%)'] = ((percentage_calculations_808['Drug Control'] - percentage_calculations_808['PBS Control 808'])
                                          / factor_control) * 100
        graph_perc_808['Drug Control'] = df['Drug Control - 808nm (%)']
        df['STD (%)'] = (percentage_calculations_808['Drug Control STD'] / factor_control) * 100
        nm_808.append(df)
    else:
        df = pd.DataFrame()
        df['Drug Control - 248nm (%)'] = ((percentage_calculations_248['Drug Control'] - 2645) / factor_control) * 100
        graph_perc_248['Drug Control'] = df['Drug Control - 248nm (%)']
        df['STD (%)'] = (percentage_calculations_248['Drug Control STD'] / factor_control) * 100
        nm_248.append(df)

        df = pd.DataFrame()
        df['Drug Control - 808nm (%)'] = ((percentage_calculations_808['Drug Control'] - 2645) / factor_control) * 100
        graph_perc_808['Drug Control'] = df['Drug Control - 808nm (%)']
        df['STD (%)'] = (percentage_calculations_808['Drug Control STD'] / factor_control) * 100
        nm_808.append(df)

    values_248 = pd.concat(nm_248, axis=1)
    values_808 = pd.concat(nm_808, axis=1)

    graph_nm_248, graph_perc_248 = convert_to_numeric_dataframes(False,
                                                         graph_nm_248, graph_perc_248)

    graph_perc_808, graph_perc_808 = convert_to_numeric_dataframes(False,
                                                                 graph_perc_808, graph_perc_808)

    plt.figure()
    for column in graph_nm_248.columns:
        if column != 'Time':
            plt.scatter(graph_nm_248['Time'], graph_nm_248[column], label=column)
            z = np.polyfit(graph_nm_248['Time'], graph_nm_248[column], 4)
            p = np.poly1d(z)
            plt.plot(sorted(graph_nm_248['Time']), p(sorted(graph_nm_248['Time'])), linestyle='--')

    plt.title(f'{drug_name} Release (nM)')
    plt.xlabel('Time (h)')
    plt.ylabel('Average (nM)')
    plt.legend()
    plt.grid(True)

    plt.savefig(f'{drug_name} Release (nM) - 248nm.png', bbox_inches='tight')
    plt.close()

    plt.figure()

    for column in graph_perc_248.columns:
        if column != 'Time':
            plt.scatter(graph_nm_248['Time'], graph_nm_248[column], label=column)
            z = np.polyfit(graph_nm_248['Time'], graph_nm_248[column], 4)
            p = np.poly1d(z)
            plt.plot(sorted(graph_nm_248['Time']), p(sorted(graph_nm_248['Time'])), linestyle='--')

    plt.title(f"{drug_name} Release (%)")
    plt.xlabel('Time (h)')
    plt.ylabel('Average (%)')
    plt.legend()
    plt.grid(True)

    plt.savefig(f'{drug_name} Release (%) - 248nm.png', bbox_inches='tight')
    plt.close()

    plt.figure()
    for column in graph_nm_808.columns:
        if column != 'Time':
            plt.scatter(graph_nm_808['Time'], graph_nm_808[column], label=column)
            z = np.polyfit(graph_nm_808['Time'], graph_nm_808[column], 4)
            p = np.poly1d(z)
            plt.plot(sorted(graph_nm_808['Time']), p(sorted(graph_nm_808['Time'])), linestyle='--')

    plt.title(f'{drug_name} Release (nM)')
    plt.xlabel('Time (h)')
    plt.ylabel('Average (nM)')
    plt.legend()
    plt.grid(True)

    plt.savefig(f'{drug_name} Release (nM) - 808nm.png', bbox_inches='tight')
    plt.close()

    plt.figure()

    for column in graph_perc_808.columns:
        if column != "Time":
            plt.scatter(graph_perc_808['Time'], graph_perc_808[column], label=column)
            z = np.polyfit(graph_perc_808['Time'], graph_perc_808[column], 4)
            p = np.poly1d(z)
            plt.plot(sorted(graph_perc_808['Time']), p(sorted(graph_perc_808['Time'])), linestyle='--')

    plt.title(f'{drug_name} Release (%)')
    plt.xlabel('Time (h)')
    plt.ylabel('Average (%)')
    plt.legend()
    plt.grid(True)

    plt.savefig(f'{drug_name} Release (%) - 808nm.png', bbox_inches='tight')
    plt.close()


    return values_248, values_808


if __name__ == "__main__":
    required_packages = ['pandas', 'openpyxl', 'xlsxwriter', 'numpy', 'matplotlib', 'sympy']
    for package in required_packages:
        if not is_module_installed(package):
            print(f"The package '{package}' is not installed. Please install it and try again.")
            exit(1)
    if len(sys.argv) < 3:
        print("Usage: python resiquimod_calculator.py <excel sheet location(str)> <sheet_name(str)>")
        sys.exit(1)
    excel_sheet_location = sys.argv[1]
    required_sheet_name = sys.argv[2]
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
    if not isinstance(replication, int) or replication == 0:
        print("The replication value is 0 or invalid. Please check the Excel sheet and try again.")
        sys.exit(1)

    drug_name = data.iloc[0, 1]
    if not isinstance(drug_name, str):
        print("Invalid value found for drug name. Please check the Excel sheet and try again.")
        sys.exit(1)

    calculated_248, calculated_808 = calculate(data, replication, drug_name)

    pH_value = data.iloc[5, 1]
    try:
        if pH_value == 0:
            print("The pH_value value is 0 or invalid. Please check the Excel sheet and try again.")
            sys.exit(1)
    except TypeError:
        print("The pH_value value is 0 or invalid. Please check the Excel sheet and try again.")
        sys.exit(1)


    output_file = f"./output/{drug_name} Release pH={pH_value}.xlsx"

    if not os.path.exists('./output'):
        os.makedirs('./output')

    if os.path.exists(output_file):
        os.remove(output_file)

    # Now write the DataFrames to an Excel file
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        # Write each DataFrame to a different sheet or the same sheet in Excel
        calculated_248.to_excel(writer, sheet_name="Absorption 248nm", startrow=0, startcol=0, index=False)
        calculated_808.to_excel(writer, sheet_name="Absorption 808nm", startrow=0, startcol=0, index=False)

        # Access the xlsxwriter workbook and worksheet objects
        workbook = writer.book
        worksheet_248 = writer.sheets["Absorption 248nm"]
        worksheet_808 = writer.sheets["Absorption 808nm"]

        # Calculate the position where you want to insert your plots
        # For example, below the DataFrames
        image_row_248 = len(calculated_248) + 4
        image_row_808 = len(calculated_808) + 4

        # Insert the plots
        worksheet_248.insert_image(image_row_248, 0, f'{drug_name} Release (nM) - 248nm.png')
        worksheet_248.insert_image(image_row_248, 4, f'{drug_name} Release (%) - 248nm.png')

        worksheet_808.insert_image(image_row_808, 0, f'{drug_name} Release (nM) - 808nm.png')
        worksheet_808.insert_image(image_row_808, 4, f'{drug_name} Release (%) - 808nm.png')


        for column in range(calculated_248.shape[1]):
            worksheet_248.set_column(column, column, 22)
        for column in range(calculated_808.shape[1]):
            worksheet_808.set_column(column, column, 22)

    os.remove(f'{drug_name} Release (nM) - 248nm.png')
    os.remove(f'{drug_name} Release (%) - 248nm.png')
    os.remove(f'{drug_name} Release (nM) - 808nm.png')
    os.remove(f'{drug_name} Release (%) - 808nm.png')

    print(f"Data has been written to {output_file}")
    print("The program will now exit.")
    sys.exit(0)
