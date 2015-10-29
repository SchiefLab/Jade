import pandas




class GeneralPandasDataFrame(pandas.DataFrame):
    def __init__(self, data=None, index=None, columns=None, dtype=None,copy=False):
        pandas.DataFrame.__init__(data=data, index=index, columns=columns, dtype=dtype,copy=copy)


