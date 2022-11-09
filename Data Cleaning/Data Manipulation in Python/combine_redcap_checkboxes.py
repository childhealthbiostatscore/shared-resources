"""
This function takes as its arguments a pandas dataframe and a dictionary linking the current column names to new names. This assumes standard REDCap formatting where the base variable is split into multiple columns, named according to the format: "base_name___*" where * is a digit (1 through the number of checkbox options). The list of levels names must be in the correct order.
"""


def combine_redcap_checkboxes(df, base_name="", levels=[], sep=" & ",
                              drop=True):
    # Libraries
    import pandas as pd
    import numpy as np
    # A list of integer strings (one for each level)
    num = [str(s) for s in range(1, len(levels) + 1)]
    # Get column names
    cols = [base_name + "___" + s for s in num]
    # Which columns are checked
    values = df[cols].apply(lambda x: " & ".join(
        [str(s + 1) for s in np.where((x == 1) | (x == "Checked"))[0]]), axis=1)
    # Dictionary for number to level translation
    column_translation = dict(zip(num, levels))
    # Replace numbers with levels
    for word, replacement in column_translation.items():
        values = [s.replace(word, replacement) for s in values]
    # Add to dataframe
    df[base_name] = values
    # Drop old columns
    if drop:
        df.drop(cols, axis=1, inplace=True)
    # Return df
    return df
