import pandas




class GeneralPandasDataFrame(pandas.DataFrame):
    def __init__(self, data=None, index=None, columns=None, dtype=None,copy=False):
        pandas.DataFrame.__init__(data=data, index=index, columns=columns, dtype=dtype,copy=copy)

    def drop_duplicate_columns(self):
        """
        Drop Duplicate columns from the DataFrame in place
        :return:
        """
        self = self.T.groupby(level=0).first().T


def drop_duplicate_columns(df):
    """
    Drop Duplicate columns from the DataFrame

    :param df: pandas.DataFrame
    :rtype: pandas.DataFrame
    """
    return df.T.groupby(level=0).first().T

