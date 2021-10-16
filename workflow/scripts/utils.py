#!/usr/bin/env python3
import pandas as pd
import sys
import h5py
import os

def make_better(name):
    # excel does not allow some special characters in sheet name
    for c in '[]:*?/\\':
        name = name.replace(c, '-')
    return name


class DF2Excel:
    def __init__(self, excel_path):
        self.writer = pd.ExcelWriter(excel_path, engine='xlsxwriter')

    def __del__(self):
        self.writer.save()

    def close(self):
        self.writer.save()

    def single_idx_df_to_excel(self, dataframe, sheet_name, write_header,
                               write_index, start_row = 0):
        # make better sheet name for excel
        sheet_name = make_better(sheet_name)
        # now write the actual data without the header
        dataframe.to_excel(
            self.writer,
            sheet_name=sheet_name,
            startrow=start_row,
            header=write_header,
            index=write_index
        )
        # create excel auto filter
        if write_header:
            self.writer.sheets[sheet_name].autofilter(
                start_row, 0, start_row, dataframe.shape[1]
            )

    def multi_idx_df_to_excel(self, dataframe, sheet_name, start_row = 0):
        # make better sheet name for excel
        sheet_name = make_better(sheet_name)
        # first write the header structure, keep bare-minimum 1 row for index
        dataframe.iloc[0:1,:].to_excel(self.writer, sheet_name=sheet_name,
                                       startrow=start_row)

        # now write the actual data without the header
        start_row +=  dataframe.columns.nlevels + 1
        freeze_row = start_row
        freeze_col = dataframe.index.nlevels
        dataframe.to_excel(
            self.writer,
            sheet_name=sheet_name,
            startrow=start_row,
            header=False,
            merge_cells=False,
            freeze_panes=(freeze_row, freeze_col)
        )
        # create excel auto filter
        self.writer.sheets[sheet_name].autofilter(
            freeze_row-1, 0,
            freeze_row-1 + dataframe.shape[0],
            freeze_col-1 + dataframe.shape[1]
        )
        # make column index names right justified
        cell_format = self.writer.book.add_format()
        cell_format.set_align('right')
        cell_format.set_align('vcenter')
        cell_format.set_bold()
        for row in range(dataframe.columns.nlevels):
            self.writer.sheets[sheet_name].write(
                row, freeze_col-1, dataframe.columns.names[row], cell_format
            )

    def write(self, dataframe, sheet_name, write_header=True,
              write_index=False, start_row=0):
        if isinstance(dataframe.columns, pd.MultiIndex):
            self.multi_idx_df_to_excel(dataframe, sheet_name, start_row)
        else:
            self.single_idx_df_to_excel(dataframe, sheet_name, write_header,
                                        write_index, start_row)

    def convert_single(self, hdf_path, key):
        df = pd.read_hdf(hdf_path, key=key)
        self.write(df, key)

    def convert(self, hdf_path):
        for key in h5py.File(hdf_path, 'r').keys():
            self.convert_single(hdf_path, key)


def pvalstr(x):
    if x < 0.001: return '***'
    if x < 0.01: return '**'
    if x < 0.05: return '*'
    return ''

def print_pass(df):
    print(df)
    print(df.columns)
    return df

def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if directory and not os.path.exists(directory):
        print("Creating new directory", directory)
        os.makedirs(directory)
    return file_path
