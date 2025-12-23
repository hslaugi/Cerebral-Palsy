import pandas as pd
from scipy import stats
from itertools import combinations

df = pd.read_csv('Results_TD_CP.csv')

metrics = ['dti_fa', 'md', 'ad', 'rd']
groups = ['TD_ND', 'TD_D', 'CSC_LA', 'CSC_MA', 'PV_LA', 'PV_MA', 'Others_LA', 'Others_MA']

# Defining Strict Threshold (Bonferroni Correction)
# We have 8 groups, so there are 28 unique pairs.
# 0.05 / 28 = 0.00178. P-values must be lower than this to be significant.
# This prevents "false positives" from running so many tests.

pairs = list(combinations(groups, 2))
corrected_alpha = 0.05 / len(pairs)

print(f"Post-Hoc Analysis (Bonferroni Alpha: {corrected_alpha:.5f})")
print("-" * 40)

for metric in metrics:
    print(f"\n {metric} Significant Pairs")
    subset = df[df['Tracts Name'] == metric]
    
    significant_found = False
    
    # Testing every possible pair of groups
    for group1, group2 in pairs:
        data1 = subset[group1].dropna().values
        data2 = subset[group2].dropna().values
        
        if len(data1) == 0 or len(data2) == 0:
            continue
            
        # Mann-Whitney U Test
        # 'two-sided' is the standard assumption
        stat, p_val = stats.mannwhitneyu(data1, data2, alternative='two-sided')
        
        if p_val < corrected_alpha:
            print(f"{group1} vs {group2} | p={p_val:.2e}")
            significant_found = True
            
    if not significant_found:
        print("No pairs met the strict significance threshold.")