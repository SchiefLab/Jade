import basic.pandas




class GeneralPandasDataFrame(basic.pandas.DataFrame):
    def __init__(self, data=None, index=None, columns=None, dtype=None,copy=False):
        basic.pandas.DataFrame.__init__(data=data, index=index, columns=columns, dtype=dtype,copy=copy)

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
