# Drug Release Calculator

This is a script which takes 3 inputs
1. The location of the excel file containing the data
2. Name of the sheet in the excel file
3. If the output should be condensed (optional, default is True)

The script will then output an excel sheet containing graphs, standard deviation
and other values required to interpret the data.

An example of the input file is given in the root folder.

## Installation
Ensure that you have python3 and pip3 installed on your system.
Once you clone this repository, run the following commands in the directory of the script to install the required packages.
1. `pip3 install -r requirements.txt`

## Usage
To run the script, use the following command
1. `python3 resiquimod_calculator.py <input_file> <sheet_name> <condense>`

Note that the condense argument is optional and will default to True if not given.

## Output
The script will output an Excel file in the directory of the script. This file will contain the graphs and data required to interpret the data.

## Example
1. `python3 resiquimod_calculator.py example.xlsx Sheet1 False`
2. `python3 resiquimod_calculator.py example.xlsx Sheet1`

If none of the above commands work, try removing the "3" from the python and pip commands.

## Footnote
Please ensure that the data is in the correct format. The script will not work if the data is not in the correct format.
I have provided an example file in the root directory to show the correct format.

This was created for the [Smith Lab](https://smithlab.iq.msu.edu/) at the [Institute of Qualitative Health Science and Engineering](https://iq.msu.edu/).