# Function for calculating "GFR_Schwartz","GFR_FAS","GFR_Zappitelli","GFR_CKIDU25" based on Laura's code for Petter Bjornstad's eGFR vs mGFR study male and female arguments specify how sex is coded.

def calc_egfr(df, age="age", serum_creatinine="creatinine_s",
              cystatin_c="cystatin_c_s", bun="bun", height="height", sex="sex",
              male="Male", female="Female", alpha=0.5):
    import pandas as pd
    import numpy as np
    # Format input
    sex = df[sex]
    sex.replace({male: "M", female: "F", "": np.nan}, inplace=True)
    qcr = np.floor(df[age])
    serum_creatinine = pd.to_numeric(df[serum_creatinine], errors="coerce")
    cystatin_c = pd.to_numeric(df[cystatin_c], errors="coerce")
    height = pd.to_numeric(df[height], errors="coerce")
    bun = pd.to_numeric(df[bun], errors="coerce")
    # Younger participants
    qcr.replace({8: 0.46, 9: 0.49, 10: 0.51, 11: 0.53,
                12: 0.57, 13: 0.59, 14: 0.61}, inplace=True)
    # Females
    qcr.loc[(qcr == 15) & (sex == "F")] = 0.64
    qcr.loc[(qcr == 16) & (sex == "F")] = 0.67
    qcr.loc[(qcr == 17) & (sex == "F")] = 0.69
    qcr.loc[(qcr == 18) & (sex == "F")] = 0.69
    qcr.loc[(qcr >= 19) & (sex == "F")] = 0.70
    # Males
    qcr.loc[(qcr == 15) & (sex == "M")] = 0.72
    qcr.loc[(qcr == 16) & (sex == "M")] = 0.78
    qcr.loc[(qcr == 17) & (sex == "M")] = 0.82
    qcr.loc[(qcr == 18) & (sex == "M")] = 0.85
    qcr.loc[(qcr == 19) & (sex == "M")] = 0.88
    qcr.loc[(qcr > 19) & (sex == "M")] = 0.90
    # Calculate final metrics
    eGFR_fas_cr = 107.3 / (serum_creatinine / qcr)
    # eGFR FAS combined creatinine and cystatin-C
    f1 = serum_creatinine / qcr
    f2 = 1 - alpha
    f3 = cystatin_c / 0.82
    eGFR_fas_cr_cysc = 107.3 / ((0.5 * f1) + (f2 * f3))
    # eGFR Zapatelli
    eGFR_Zap = (507.76 * np.exp((0.003 * height)) /
                ((cystatin_c ** 0.635) * ((serum_creatinine * 88.4) ** 0.547)))
    # eGFR Schwartz
    m = sex.replace({"M": 1, "F": 0})
    eGFR_Schwartz = 39.1 * ((height / serum_creatinine) ** 0.516) * ((1.8 / cystatin_c) ** 0.294) * ((30 / bun) ** 0.169) * \
        (1.099 ** m) * ((height / 1.4) ** 0.188)
    # eGFR bedside Schwartz
    eGFR_bedside_Schwartz = (41.3 * (height / 100)) / serum_creatinine
    # Combine
    egfr = pd.concat(
        [eGFR_Schwartz, eGFR_bedside_Schwartz, eGFR_Zap, eGFR_fas_cr, eGFR_fas_cr_cysc], axis=1)
    egfr.columns = ["eGFR_Schwartz", "eGFR_bedside_Schwartz", "eGFR_Zap",
                    "eGFR_fas_cr", "eGFR_fas_cr_cysc"]
    df = pd.concat([df, egfr], axis=1)
    # Return
    return df
