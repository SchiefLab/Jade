import pandas as pd
import os



class GeneralPandasDataFrame(pd.DataFrame):
    def __init__(self, data=None, index=None, columns=None, dtype=None,copy=False):
        pd.DataFrame.__init__(data=data, index=index, columns=columns, dtype=dtype,copy=copy)

    def drop_duplicate_columns(self):
        """
        Drop Duplicate columns from the DataFrame in place
        :return:
        """

        #I'm not sure how to do this inplace, without reassigning self.  If you know, please edit this.


        self = self.T.groupby(level=0).first().T

    def detect_numeric(self):
        self = self.convert_objects(convert_numeric=True)

    def get_columns(self, columns):
        return self[columns]

    def get_matches(self, column, to_match):
        """
        Get all the rows that match a paricular element of a column.
        :param column: str
        :param to_match: str
        :rtype: pandas.DataFrame
        """

        return self[self[column] == to_match]

    def get_row_matches(self, column1, to_match, column2):
        """
        Get the elements of the rows that match a particular column.  If one element, this can be converted easily enough
        :param column1: str
        :param to_match: str
        :param column2: str
        :rtype: pandas.Series
        """

        return self[self[column1] == to_match][column2]

def multi_tab_excel(df_list, sheet_list, file_name):
    """
    Writes multiple dataframes as separate sheets in an output excel file.

    If directory of output does not exist, it will create it.

    Author: Tom Dobbs
    http://stackoverflow.com/questions/32957441/putting-many-python-pandas-dataframes-to-one-excel-worksheet


    :param df_list: [pd.Dataframe]
    :param sheet_list: [str]
    :param file_name: str

    """
    if not os.path.exists(os.path.dirname(file_name)):
        os.mkdir(os.path.dirname(file_name))

    writer = pd.ExcelWriter(file_name,engine='xlsxwriter')
    for dataframe, sheet in zip(df_list, sheet_list):
        dataframe.to_excel(writer, sheet_name=sheet, startrow=0 , startcol=0)
    writer.save()


def drop_duplicate_columns(df):
    """
    Drop Duplicate columns from the DataFrame.
    Return DF

    :param df: pandas.DataFrame
    :rtype: pandas.DataFrame
    """
    return df.T.groupby(level=0).first().T

def detect_numeric(df):
    """
    Detect numeric components

    :param df: pandas.DataFrame
    :rtype: pandas.DataFrame

    """
    return df.convert_objects(convert_numeric=True)

def get_columns(df, columns):
    """
    Get a new dataframe of only the columns
    :param df: pandas.DataFrame
    :param columns: list
    :rtype: pandas.DataFrame
    """
    return df[columns]

def get_matches(df, column, to_match):
    """
    Get all the rows that match a paricular element of a column.
    :param df: pandas.DataFrame
    :param column: str
    :param to_match: str
    :rtype: pandas.DataFrame
    """

    return df[df[column] == to_match]

def get_match_by_array(df, column, match_array):
    """
    Get a new dataframe of all dataframes of the subset series, match_array

    Note: This will result in a dataframe, but there may be strange issues when you go to plot the data in seaborn
            No idea why.

    :param df: pandas.DataFrame
    :param column: str
    :param match_array: pandas.Series
    :rtype: pandas.DataFrame
    """

    new_df = df[df[column].isin(match_array)]
    return new_df


def get_row_matches(df, column1, to_match, column2):
    """
    Get the elements of the rows that match a particular column.  If one element, this can be converted easily enough
    :param df: pandas.DataFrame
    :param column1: str
    :param to_match: str
    :param column2: str
    :rtype: pandas.Series
    """

    return df[df[column1] == to_match][column2]
