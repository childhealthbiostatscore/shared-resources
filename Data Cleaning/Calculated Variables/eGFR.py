# Function for calculating "GFR_Schwartz","GFR_FAS","GFR_Zappitelli","GFR_CKIDU25" based on Laura's code for Petter Bjornstad's eGFR vs mGFR study male and female arguments specify how sex is coded.

def calc_egfr(df, age="age", serum_creatinine="creatinine_s",
              cystatin_c="cystatin_c_s", bun="bun", height="height", sex="sex",
              male="Male", female="Female", alpha=0.5):
    import pandas as pd
    import numpy as np
    # Format input
    df[sex].replace({male: "M", female: "F", "": np.nan}, inplace=True)
    df["age_floor"] = np.floor(df[age])
    df["qcr"] = np.nan
    df[serum_creatinine] = pd.to_numeric(df[serum_creatinine], errors="coerce")
    df[cystatin_c] = pd.to_numeric(df[cystatin_c], errors="coerce")
    df[height] = pd.to_numeric(df[height], errors="coerce")
    df[bun] = pd.to_numeric(df[bun], errors="coerce")
    # Younger participants
    df.loc[df["age_floor"] == 8, "qcr"] = 0.46
    df.loc[df["age_floor"] == 9, "qcr"] = 0.49
    df.loc[df["age_floor"] == 10, "qcr"] = 0.51
    df.loc[df["age_floor"] == 11, "qcr"] = 0.53
    df.loc[df["age_floor"] == 12, "qcr"] = 0.57
    df.loc[df["age_floor"] == 13, "qcr"] = 0.59
    df.loc[df["age_floor"] == 14, "qcr"] = 0.61
    # Females
    df.loc[(df["age_floor"] == 15) & (df[sex] == "F"), "qcr"] = 0.64
    df.loc[(df["age_floor"] == 16) & (df[sex] == "F"), "qcr"] = 0.67
    df.loc[(df["age_floor"] == 17) & (df[sex] == "F"), "qcr"] = 0.69
    df.loc[(df["age_floor"] == 18) & (df[sex] == "F"), "qcr"] = 0.69
    df.loc[(df["age_floor"] == 19) & (df[sex] == "F"), "qcr"] = 0.70
    df.loc[(df["age_floor"] > 19) & (df[sex] == "F"), "qcr"] = 0.70
    # Males
    df.loc[(df["age_floor"] == 15) & (df[sex] == "M"), "qcr"] = 0.72
    df.loc[(df["age_floor"] == 16) & (df[sex] == "M"), "qcr"] = 0.78
    df.loc[(df["age_floor"] == 17) & (df[sex] == "M"), "qcr"] = 0.82
    df.loc[(df["age_floor"] == 18) & (df[sex] == "M"), "qcr"] = 0.85
    df.loc[(df["age_floor"] == 19) & (df[sex] == "M"), "qcr"] = 0.88
    df.loc[(df["age_floor"] > 19) & (df[sex] == "M"), "qcr"] = 0.90
    # Calculate final metrics
    df["eGFR_fas_cr"] = 107.3 / (df[serum_creatinine] / df["qcr"])
    # eGFR FAS combined creatinine and cystatin-C
    f1 = df[serum_creatinine] / df["qcr"]
    f2 = 1 - alpha
    f3 = df[cystatin_c] / 0.82
    df["eGFR_fas_cr_cysc"] = 107.3 / ((0.5 * f1) + (f2 * f3))
    # eGFR Zapatelli
    df["eGFR_Zap"] = (507.76 * np.exp(0.003 * (df[height]))) / \
        ((df[cystatin_c] ** 0.635) * ((df[serum_creatinine] * 88.4) ** 0.547))
    # eGFR Schwartz
    df["m"] = df[sex].replace({"M": 1, "F": 0})
    df["eGFR_Schwartz"] = 39.1 * ((df[height] / df[serum_creatinine]) ** 0.516) * ((1.8 / df[cystatin_c]) ** 0.294) * ((30 / df[bun]) ** 0.169) * \
        (1.099 ** df["m"]) * ((df[height] / 1.4) ** 0.188)
    # eGFR bedside Schwartz
    df["eGFR_bedside_Schwartz"] = (
        41.3 * (df[height] / 100)) / df[serum_creatinine]
    # remove calculation columns
    df.drop(["m", "age_floor"], axis=1, inplace=True)
    # Return
    return df
