import os
import sys

import pandas as pd
import openpyxl, xlsxwriter, re, pyarrow
import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, Eq, solve
from enum import Enum
import importlib.util

class Type(Enum):
    FIRST_PEAK = 1
    SECOND_PEAK = 2

# Check if the variable is of the desired type
def variable_check(var, desired_type):
    if not isinstance(var, desired_type):
        print(f"Invalid value found for {var}. Please check the Excel sheet and try again.")
        sys.exit(1)

# Read the Excel sheet
def read_excel_sheet(location) -> (pd.ExcelFile, dict):
    # Read the Excel sheet and return the ExcelFile and the sheets
    xl = pd.ExcelFile(location)

    # Get the sheet names
    sheet_names = xl.sheet_names

    sheets = {}

    # Read the sheets
    for sheet in sheet_names:
        sheets[sheet] = xl.parse(sheet, header=None)

    # Return the ExcelFile and the sheets
    return xl, sheets

# Check if the module is installed
def is_module_installed(module_name) -> bool:
    spec = importlib.util.find_spec(module_name)
    return spec is not None

# Calculate the standard deviation
def calculate_std(row, nm_type, *args):
    # Apply the operation to specific cells and calculate the standard deviation
    if nm_type == Type.FIRST_PEAK:
        return np.std([((cell + args[0]) / args[1]) for cell in row], ddof=0)
    elif nm_type == Type.SECOND_PEAK:
        return np.std([(cell * 1000000) / 7900 for cell in row], ddof=0)


# Get the slope and intercept
def get_slope_intercept(equation):
    equation = equation.replace("^", "**")
    equation = equation.replace("y=", "")
    x, y = symbols('x y')  # Define 'x' and 'y'
    eq = Eq(y, eval(equation))
    solution = solve(eq, y)
    slope = solution[0].coeff(x)
    intercept = solution[0].subs(x, 0)
    intercept = abs(intercept)
    return float(slope), float(intercept)

# Reset the indexes and convert the data to numeric
def reset_and_convert(exit_if_null, *data) -> tuple[pd.DataFrame, ...]:
    result = []
    for i in range(len(data)):
        df = data[i]
        # Check if the data is a Series
        if isinstance(df, pd.Series):
            # Convert the Series to a DataFrame
            df = df.to_frame()

        # Reset the index
        df = df.reset_index(drop=True)
        if df.shape[1] == 1:
            df.columns = [0]
        else:
            df.columns = range(df.shape[1])

        # Convert the data to numeric
        new_data = df.apply(pd.to_numeric, errors='coerce')
        if exit_if_null:
            # Check if the data has any null values
            if new_data.isnull().values.any():
                print("Please check the data provided. There are blank cells in the data.")
                print(new_data)
                print(f"The argument number is {i + 1}. Please check the Excel sheet and try again.")
                sys.exit(1)
        new_data = new_data.dropna()
        result.append(new_data)
    return tuple(result)

# Calculate the concentration factors
def concentration_calculations(data, loaded, control) -> tuple[float, float]:
    interior_volume = data.iloc[3, 1]
    exterior_volume = data.iloc[4, 1]
    loaded_swnts = (loaded * (interior_volume / exterior_volume))
    control_swnts = (control * (interior_volume / exterior_volume))
    return loaded_swnts, control_swnts

# Find the rows for the First and Second Peaks
def find_row(data, nanometer, first_peak, second_peak) -> slice:
    start_row, end_row = 0, 0
    time_column = data.iloc[:, 3]
    if nanometer == Type.FIRST_PEAK:
        for index, value in enumerate(time_column):
            if value == f'Absorption {first_peak} nm':
                start_row = index + 2
            if value == f'Absorption {second_peak} nm':
                end_row = index - 1
                break
        return slice(start_row, end_row)
    elif nanometer == Type.SECOND_PEAK:
        for index, value in enumerate(time_column):
            if value == f'Absorption {second_peak} nm':
                start_row = index + 2
        end_row = time_column.size
        return slice(start_row, end_row)

# Find the columns for the data
def find_columns(data, drug_name) -> dict[str, list]:
    header_row = data.iloc[1, :17]
    cols = {"Time": [], drug_name: [], "SWNTs Control": [], "Drug Control": [], "DDI Water Control": []}
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
        elif value.startswith("DDI Water"):
            cols["DDI Water Control"].append(index)
    return cols

# Find the peaks
def find_peaks(data) -> tuple[str, str]:
    first_peak = data.iloc[0, 3]
    match = re.search(r'Absorption (\d+) nm', first_peak)
    number = match.group(1) if match else None
    peak_column = data.iloc[1:, 3]
    try:
        second_peak = peak_column.to_string()
    except Exception as error:
        print(f"Invalid data provided. {error}. Check the Excel sheet and try again.")
        sys.exit(1)

    # Check if the second_peak has been found
    if not second_peak:
        print("Invalid data provided. No second peak found. Check the Excel sheet and try again.")
        sys.exit(1)

    match = re.search(r'Absorption (\d+) nm', second_peak)
    number2 = match.group(1) if match else None

    if number is None or number2 is None:
        print("Invalid data provided. No peaks found. Check the Excel sheet and try again.")
        sys.exit(1)

    return number, number2

# Graph the data
def graph_data(time, graph_data, title, x_label, y_label, file_name) -> None:
    plt.figure()
    for column in graph_data.columns:
        if column.endswith("STD"):
            continue
        std_column = column + " - STD"
        if std_column not in graph_data.columns:
            print(f"Warning: No STD column found for {column}")
            continue
        plt.errorbar(time, graph_data[column], yerr=graph_data[std_column], label=column, fmt='o', linestyle='dashed')
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend()
    plt.grid(True)

    plt.savefig(file_name, bbox_inches='tight')
    plt.close()

def calculate(data, replication, drug_name, first_peak, second_peak) -> tuple[pd.DataFrame, pd.DataFrame]:
    # Create Dataframe List for each peak
    first_peak_dfs = []
    second_peak_dfs = []

    # Calculate the factor for the loaded SWNTs and the control SWNTs
    factor_loaded_swnts, factor_control = concentration_calculations(data, data.iloc[9, 1], data.iloc[11, 1])

    # Calculate the slope and intercept for the equation
    slope, intercept = get_slope_intercept(data.iloc[6, 1])

    # Initialize the DataFrames for the DDI Water Control
    pd_ddi_fp, pd_ddi_sp, absorption_ddi_fp, absorption_ddi_sp, ddi_pres = None, None, None, None, False

    # Find the columns for the data
    cols = find_columns(data, drug_name)

    # Find the rows for the First Peak Data
    rows_fp = find_row(data, Type.FIRST_PEAK, first_peak, second_peak)
    if rows_fp.start == 0 and rows_fp.stop == 0:
        print(f"Invalid data provided. Invalid values found for {first_peak} nm. Check the Excel sheet and try again.")
        sys.exit(1)

    # Find the rows for the Second Peak Data
    rows_sp = find_row(data, Type.SECOND_PEAK, first_peak, second_peak)
    if rows_sp.start == 0 and rows_sp.stop == 0:
        print(f"Invalid data provided. Invalid values found for {second_peak} nm. Check the Excel sheet and try again.")
        sys.exit(1)

    # Check if the DDI Water Control data is present
    if len(cols['DDI Water Control']) != 0:
        ddi_pres = True

    # Reset and convert DDI Water Control Data
    if ddi_pres:
        absorption_ddi_fp, absorption_ddi_sp = reset_and_convert(True,
                                                                       data.iloc[rows_fp.start:rows_fp.stop, cols['DDI Water Control']],
                                                                       data.iloc[rows_sp.start:rows_sp.stop, cols['DDI Water Control']])
    # Initialize the empty DataFrames for the First and Second Peaks
    graph_nm_fp = pd.DataFrame()
    graph_perc_fp = pd.DataFrame()

    graph_nm_sp = pd.DataFrame()
    graph_perc_sp = pd.DataFrame()

    # Reset and convert the data for the First and Second Peaks
    (absorption_sp_loaded_drug, absorption_sp_swnt_control, absorption_sp_drug_control,
     absorption_fp_loaded_drug, absorption_fp_swnt_control, absorption_fp_drug_control,
     time) = reset_and_convert(True,
                               data.iloc[rows_sp.start:rows_sp.stop, cols[drug_name]],
                      data.iloc[rows_sp.start:rows_sp.stop, cols['SWNTs Control']],
                      data.iloc[rows_sp.start:rows_sp.stop, cols['Drug Control']],
                      data.iloc[rows_fp.start:rows_fp.stop, cols[drug_name]],
                      data.iloc[rows_fp.start:rows_fp.stop, cols['SWNTs Control']],
                      data.iloc[rows_fp.start:rows_fp.stop, cols['Drug Control']],
                      data.iloc[rows_fp.start:rows_fp.stop, cols['Time']])

    # Add the Time column to the DataFrames
    time.columns = ['Time (h)']
    first_peak_dfs.append(time)
    second_peak_dfs.append(time)

    # First Concentration Peak Calculation

    # Calculate Drug Concentration for the First Peak
    pd_loaded_fp = pd.DataFrame()
    pd_loaded_fp[f'{drug_name} SWNTs - {first_peak}nm (nM)'] = ((absorption_fp_loaded_drug.sum(axis=1) / replication) + intercept) / slope
    graph_nm_fp[f'{drug_name} SWNTs'] = pd_loaded_fp[f'{drug_name} SWNTs - {first_peak}nm (nM)']
    pd_loaded_fp['STD (nM)'] = absorption_fp_loaded_drug.apply(calculate_std, args=(Type.FIRST_PEAK, intercept, slope), axis=1)
    graph_nm_fp[f'{drug_name} SWNTs - STD'] = pd_loaded_fp['STD (nM)']
    first_peak_dfs.append(pd_loaded_fp)

    # Calculate SWNTs Control Concentration for the First Peak
    pd_swnt_fp = pd.DataFrame()
    pd_swnt_fp[f'SWNTs Control - {first_peak}nm (nM)'] = ((absorption_fp_swnt_control.sum(axis=1) / replication) + intercept) / slope
    graph_nm_fp['SWNTs Control'] = pd_swnt_fp[f'SWNTs Control - {first_peak}nm (nM)']
    pd_swnt_fp['STD (nM)'] = absorption_fp_swnt_control.apply(calculate_std, args=(Type.FIRST_PEAK, intercept, slope), axis=1)
    graph_nm_fp['SWNTs Control - STD'] = pd_swnt_fp['STD (nM)']
    first_peak_dfs.append(pd_swnt_fp)

    # Calculate Drug Control Concentration for the First Peak
    pd_drug_fp = pd.DataFrame()
    pd_drug_fp[f'Drug Control - {first_peak}nm (nM)'] = ((absorption_fp_drug_control.sum(axis=1) / replication) + intercept) / slope
    graph_nm_fp['Drug Control'] = pd_drug_fp[f'Drug Control - {first_peak}nm (nM)']
    pd_drug_fp['STD (nM)'] = absorption_fp_drug_control.apply(calculate_std, args=(Type.FIRST_PEAK, intercept, slope), axis=1)
    graph_nm_fp['Drug Control - STD'] = pd_drug_fp['STD (nM)']
    first_peak_dfs.append(pd_drug_fp)

    # End First Concentration Peak Calculation
    # Start Second Concentration Peak Calculation

    # Calculate Drug Concentration for the Second Peak
    pd_loaded_sp = pd.DataFrame()
    pd_loaded_sp[f'{drug_name} SWNTs - {second_peak}nm (nM)'] = ((absorption_sp_loaded_drug.sum(axis=1) / replication) * 1000000) / 7900
    graph_nm_sp[f'{drug_name} SWNTs'] = pd_loaded_sp[f'{drug_name} SWNTs - {second_peak}nm (nM)']
    pd_loaded_sp['STD (nM)'] = absorption_sp_loaded_drug.apply(calculate_std, args=(Type.SECOND_PEAK, None, None), axis=1)
    graph_nm_sp[f'{drug_name} SWNTs - STD'] = pd_loaded_sp['STD (nM)']
    second_peak_dfs.append(pd_loaded_sp)

    # Calculate SWNTs Control Concentration for the Second Peak
    pd_swnt_sp = pd.DataFrame()
    pd_swnt_sp[f'SWNTs Control - {second_peak}nm (nM)'] = ((absorption_sp_swnt_control.sum(axis=1) / replication) * 1000000) / 7900
    graph_nm_sp['SWNTs Control'] = pd_swnt_sp[f'SWNTs Control - {second_peak}nm (nM)']
    pd_swnt_sp['STD (nM)'] = absorption_sp_swnt_control.apply(calculate_std, args=(Type.SECOND_PEAK, None, None), axis=1)
    graph_nm_sp['SWNTs Control - STD'] = pd_swnt_sp['STD (nM)']
    second_peak_dfs.append(pd_swnt_sp)

    # Calculate Drug Control Concentration for the Second Peak
    pd_drug_sp = pd.DataFrame()
    pd_drug_sp[f'Drug Control - {second_peak}nm (nM)'] = ((absorption_sp_drug_control.sum(axis=1) / replication) * 1000000) / 7900
    graph_nm_sp['Drug Control'] = pd_drug_sp[f'Drug Control - {second_peak}nm (nM)']
    pd_drug_sp['STD (nM)'] = absorption_sp_drug_control.apply(calculate_std, args=(Type.SECOND_PEAK, None, None), axis=1)
    graph_nm_sp['Drug Control - STD'] = pd_drug_sp['STD (nM)']
    second_peak_dfs.append(pd_drug_sp)

    # End Second Concentration Peak Calculation
    # Start DDI Water Control Calculation

    if ddi_pres:
        # Only if DDI Water Control data is present
        # Calculate DDI Water Control Concentration for the First Peak
        pd_ddi_fp = pd.DataFrame()
        pd_ddi_fp[f'DDI Water Control - {first_peak}nm (nM)'] = ((absorption_ddi_fp.sum(axis=1) / replication) + intercept) / slope
        graph_nm_fp['DDI Water Control'] = pd_ddi_fp[f'DDI Water Control - {first_peak}nm (nM)']
        pd_ddi_fp['STD (nM)'] = absorption_ddi_fp.apply(calculate_std, args=(Type.FIRST_PEAK, intercept, slope), axis=1)
        graph_nm_fp['DDI Water Control - STD'] = pd_ddi_fp['STD (nM)']
        first_peak_dfs.append(pd_ddi_fp)

        # Calculate DDI Water Control Concentration for the Second Peak
        pd_ddi_sp = pd.DataFrame()
        pd_ddi_sp[f'DDI Water Control - {second_peak}nm (nM)'] = ((absorption_ddi_sp.sum(axis=1) / replication) * 1000000) / 7900
        graph_nm_sp['DDI Water Control'] = pd_ddi_sp[f'DDI Water Control - {second_peak}nm (nM)']
        pd_ddi_sp['STD (nM)'] = absorption_ddi_sp.apply(calculate_std, args=(Type.SECOND_PEAK, None, None), axis=1)
        graph_nm_sp['DDI Water Control - STD'] = pd_ddi_sp['STD (nM)']
        second_peak_dfs.append(pd_ddi_sp)

    # End DDI Water Control Calculation
    # Start Percentage Calculations

    # Start Non Normalized Data
    # Calculate the percentage for the First Peak
    df = pd.DataFrame()
    df[f'{drug_name} SWNTs - {first_peak} (%)'] = ((pd_loaded_fp[f'{drug_name} SWNTs - {first_peak}nm (nM)']
                                                    / factor_loaded_swnts) * 100)
    df['STD (%)'] = (pd_loaded_fp['STD (nM)'] / factor_loaded_swnts) * 100
    first_peak_dfs.append(df)

    # Calculate the percentage for the Second Peak
    df = pd.DataFrame()
    df[f'{drug_name} SWNTs - {second_peak} (%)'] = ((pd_loaded_sp[f'{drug_name} SWNTs - {second_peak}nm (nM)']
                                                     / factor_loaded_swnts) * 100)
    df['STD (%)'] = (pd_loaded_sp['STD (nM)'] / factor_loaded_swnts) * 100
    second_peak_dfs.append(df)

    # End Non Normalized Data
    # Start Normalized Data

    # Calculate the percentage for the First Peak
    df = pd.DataFrame()
    df[f'{drug_name} SWNTs - {first_peak} Normalized (%)'] = (((pd_loaded_fp[f'{drug_name} SWNTs - {first_peak}nm (nM)'] -
                               pd_swnt_fp[f'SWNTs Control - {first_peak}nm (nM)']) / factor_loaded_swnts) * 100)
    graph_perc_fp[f'{drug_name} SWNTs'] = df[f'{drug_name} SWNTs - {first_peak} Normalized (%)']
    df['STD (%)'] = (pd_loaded_fp['STD (nM)'] / factor_loaded_swnts) * 100
    graph_perc_fp[f'{drug_name} SWNTs - STD'] = df['STD (%)']
    first_peak_dfs.append(df)

    # Calculate the percentage for the Second Peak
    df = pd.DataFrame()
    df[f'{drug_name} SWNTs - {second_peak} Normalized (%)'] = (((pd_loaded_sp[f'{drug_name} SWNTs - {second_peak}nm (nM)'] -
                               pd_swnt_sp[f'SWNTs Control - {second_peak}nm (nM)']) / factor_loaded_swnts) * 100)
    graph_perc_sp[f'{drug_name} SWNTs'] = df[f'{drug_name} SWNTs - {second_peak} Normalized (%)']
    df['STD (%)'] = (pd_loaded_sp['STD (nM)'] / factor_loaded_swnts) * 100
    graph_perc_sp[f'{drug_name} SWNTs - STD'] = df['STD (%)']
    second_peak_dfs.append(df)

    # End Normalized Data

    if ddi_pres:
        # Only if DDI Water Control data is present
        # Calculate the percentage for the First Peak
        df = pd.DataFrame()
        df[f'Drug Control - {first_peak}nm Normalized (%)'] = ((pd_drug_fp[f'Drug Control - {first_peak}nm (nM)'] -
                                            pd_ddi_fp[f'DDI Water Control - {first_peak}nm (nM)'])
                                          / factor_control) * 100
        graph_perc_fp['Drug Control'] = df[f'Drug Control - {first_peak}nm Normalized (%)']
        df['STD (%)'] = (pd_drug_fp['STD (nM)'] / factor_control) * 100
        graph_perc_fp['Drug Control - STD'] = df['STD (%)']
        first_peak_dfs.append(df)

        # Calculate the percentage for the Second Peak
        df = pd.DataFrame()
        df[f'Drug Control - {second_peak}nm Normalized (%)'] = ((pd_drug_sp[f'Drug Control - {second_peak}nm (nM)']
                                                      - pd_ddi_sp[f'DDI Water Control - {second_peak}nm (nM)'])
                                          / factor_control) * 100
        graph_perc_sp['Drug Control'] = df[f'Drug Control - {second_peak}nm Normalized (%)']
        df['STD (%)'] = (pd_drug_sp['STD (nM)'] / factor_control) * 100
        graph_perc_sp['Drug Control - STD'] = df['STD (%)']
        second_peak_dfs.append(df)
    else:
        # If not, subtract 2645 from the Drug Control and divide by the factor
        # Calculate the percentage for the First Peak
        df = pd.DataFrame()
        df[f'Drug Control - {first_peak}nm (%)'] = ((pd_drug_fp[f'Drug Control - {first_peak}nm (nM)'] - 2645) / factor_control) * 100
        graph_perc_fp['Drug Control'] = df[f'Drug Control - {first_peak}nm (%)']
        df['STD (%)'] = (pd_drug_fp['STD (nM)'] / factor_control) * 100
        graph_perc_fp['Drug Control - STD'] = df['STD (%)']
        first_peak_dfs.append(df)

        # Calculate the percentage for the Second Peak
        df = pd.DataFrame()
        df[f'Drug Control - {second_peak}nm (%)'] = ((pd_drug_sp[f'Drug Control - {second_peak}nm (nM)'] - 2645) / factor_control) * 100
        graph_perc_sp['Drug Control'] = df[f'Drug Control - {second_peak}nm (%)']
        df['STD (%)'] = (pd_drug_sp['STD (nM)'] / factor_control) * 100
        graph_perc_sp['Drug Control - STD'] = df['STD (%)']
        second_peak_dfs.append(df)

    # End Percentage Calculations

    # Concatenate the DataFrames
    values_fp = pd.concat(first_peak_dfs, axis=1)
    values_sp = pd.concat(second_peak_dfs, axis=1)

    # Graph the Data
    graph_data(time,
               graph_nm_fp,
               f'{drug_name} Release (nM)',
               'Time (h)',
               'Release (nM)',
               f'{drug_name} Release (nM) - {first_peak}nm.png')

    graph_data(time,
               graph_perc_fp,
               f'{drug_name} Release (%)',
               'Time (h)',
               'Release (%)',
               f'{drug_name} Release (%) - {first_peak}nm.png')

    graph_data(time,
                graph_nm_sp,
                f'{drug_name} Release (nM)',
                'Time (h)',
                'Release (nM)',
                f'{drug_name} Release (nM) - {second_peak}nm.png')

    graph_data(time,
                graph_perc_sp,
                f'{drug_name} Release (%)',
                'Time (h)',
                'Release (%)',
                f'{drug_name} Release (%) - {second_peak}nm.png')

    # Return the DataFrames
    return values_fp, values_sp


if __name__ == "__main__":
    # Check if the required packages are installed
    required_packages = ['pandas', 'openpyxl', 'xlsxwriter', 'numpy', 'matplotlib', 'sympy']
    for package in required_packages:
        if not is_module_installed(package):
            print(f"The package '{package}' is not installed. Please install it and try again.")
            exit(1)
    # Check if the required arguments are provided
    if len(sys.argv) < 3:
        print("Usage: python drug_release_calculator.py <excel sheet location(str)> <sheet_name(str)>")
        sys.exit(1)
    excel_sheet_location = sys.argv[1]
    required_sheet_name = sys.argv[2]
    # Check if the file exists
    if not (os.path.exists(excel_sheet_location) and os.path.isfile(excel_sheet_location)):
        # Check if the file is an Excel file
        if not excel_sheet_location.endswith('.xlsx'):
            print("The file provided is not an Excel file")
            sys.exit(1)
        print("The file provided does not exist")
        sys.exit(1)
    # Read the Excel sheet
    file, full_data = read_excel_sheet(excel_sheet_location)

    # Check if the sheet name provided exists in the Excel file
    if required_sheet_name not in full_data:
        print("The sheet name provided does not exist in the Excel file")
        sys.exit(1)

    # Check if the required data is present
    data = full_data[required_sheet_name]

    # Check if the replication is valid
    replication = data.iloc[7, 1]
    if not isinstance(replication, int) or replication == 0:
        print("The replication value is 0 or invalid. Please check the Excel sheet and try again.")
        sys.exit(1)

    # Check if the drug name is valid
    drug_name = data.iloc[0, 1]
    if not isinstance(drug_name, str):
        print("Invalid value found for drug name. Please check the Excel sheet and try again.")
        sys.exit(1)

    # Check if the loaded SWNTs and control SWNTs are valid
    if data.iloc[9, 1] == 0 or data.iloc[11, 1] == 0 or variable_check(data.iloc[9, 1], int) or variable_check(data.iloc[11, 1], int):
        print("Invalid values found for loaded SWNTs and control SWNTs. Please check the Excel sheet and try again.")
        sys.exit(1)

    # Find the first and second peaks
    first_peak, second_peak = find_peaks(data)

    # Calculate the data
    calculated_fp, calculated_sp = calculate(data, replication, drug_name, first_peak, second_peak)

    # Check if the pH value is valid
    pH_value = data.iloc[5, 1]
    try:
        if pH_value == 0:
            print("The pH_value value is 0 or invalid. Please check the Excel sheet and try again.")
            sys.exit(1)
    except TypeError:
        print("The pH_value value is 0 or invalid. Please check the Excel sheet and try again.")
        sys.exit(1)


    output_file = f"./output/{drug_name} Release pH={pH_value}.xlsx"

    # Check if the output directory exists
    if not os.path.exists('./output'):
        os.makedirs('./output')

    # Check if the file exists and remove it
    if os.path.exists(output_file):
        os.remove(output_file)

    # Now write the DataFrames to an Excel file
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        # Write each DataFrame to a different sheet or the same sheet in Excel
        calculated_fp.to_excel(writer, sheet_name=f"Absorption {first_peak}nm", startrow=0, startcol=0, index=False)
        calculated_sp.to_excel(writer, sheet_name=f"Absorption {second_peak}nm", startrow=0, startcol=0, index=False)

        # Access the xlsxwriter workbook and worksheet objects
        workbook = writer.book
        worksheet_fp = writer.sheets[f"Absorption {first_peak}nm"]
        worksheet_sp = writer.sheets[f"Absorption {second_peak}nm"]

        # Calculate the position where you want to insert your plots
        # For example, below the DataFrames
        image_row_fp = len(calculated_fp) + 4
        image_row_sp = len(calculated_sp) + 4

        # Insert the plots
        worksheet_fp.insert_image(image_row_fp, 0, f'{drug_name} Release (nM) - {first_peak}nm.png')
        worksheet_fp.insert_image(image_row_fp, 4, f'{drug_name} Release (%) - {first_peak}nm.png')

        worksheet_sp.insert_image(image_row_sp, 0, f'{drug_name} Release (nM) - {second_peak}nm.png')
        worksheet_sp.insert_image(image_row_sp, 4, f'{drug_name} Release (%) - {second_peak}nm.png')

        # Set the column width and format
        for column in range(calculated_fp.shape[1]):
            worksheet_fp.set_column(column, column, 25)
        for column in range(calculated_sp.shape[1]):
            worksheet_sp.set_column(column, column, 25)

    # Remove the images
    os.remove(f'{drug_name} Release (nM) - {first_peak}nm.png')
    os.remove(f'{drug_name} Release (%) - {first_peak}nm.png')
    os.remove(f'{drug_name} Release (nM) - {second_peak}nm.png')
    os.remove(f'{drug_name} Release (%) - {second_peak}nm.png')

    # Print the output file
    print(f"Data has been written to {output_file}")
    print("The program will now exit.")
    sys.exit(0)
