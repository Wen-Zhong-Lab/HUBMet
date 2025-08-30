import os
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import pyreadr

# creadt folder for output
os.makedirs("../Output/M2_MSEA", exist_ok=True)

# load data
db = pyreadr.read_r("../Data/HUBMet_term.RDS")
HUBMet_term = db[None]
 

def table_export(da, database):
    if database == "smpdb":
        da = da.rename(columns={da.columns[0]: "SMPDB.Name", 
                                da.columns[1]: "SMPDB.ID", 
                                da.columns[2]: "Pathway.Class"})
    elif database == "kegg":
        da = da.rename(columns={da.columns[0]: "KEGG.Name", 
                                da.columns[1]: "KEGG.ID", 
                                da.columns[2]: "Pathway.Class"})
    elif database == "reactome":
        da = da.rename(columns={da.columns[0]: "Reactome.Name", 
                                da.columns[1]: "Reactome.ID"})
        da = da.drop(da.columns[[2, 3]], axis=1)
    elif database == "humanGEM":
        da = da.rename(columns={da.columns[0]: "HumanGEM-Subsystem"})
        da = da.drop(da.columns[1], axis=1)
    elif database in ["class", "subclass", "superclass", "disease"]:
        da = da.rename(columns={da.columns[0]: f"HMDB-{da.columns[0].title()}"})
    elif database in ["HBMD","HUBMet","Drug", "Disease", "Class", "Pathway"]:
        database = f"HUBMet_{database}"
    
    res_output = da.drop(columns=[col for col in ["HBM_ID", "HMDB.ID", "FDR_pri"] if col in da.columns])

    res_output = res_output.rename(columns={
        "N_met": "TermSize",
        "match_N": "N_Hit",
        "hit_n": "N_Hit",
        "match": "Hit",
        "hit_HBM_ID": "Hit",
        "hit_HMDB_ID": "Hit"
    })

    output_path = f"../Output/M2_MSEA/{database}_MSEA.txt"
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    res_output.to_csv(output_path, sep='\t', index=False)   
 


def calculate_ES_fast(args):
    # this functon returns a enrichment score of one of the 100 permetation tests
    sd, da_per, overlap_n_per, met_id_n_per = args
    np.random.seed(sd)
    
    # Generate indices for hits
    hit_indices = np.random.choice(da_per.index, overlap_n_per, replace=False)
    # hit mask
    hit_mask = np.zeros(len(da_per), dtype=bool)
    hit_mask[hit_indices] = True
    # Assign hits and calculate scores
    met_id_fc = da_per['met_id_fc'].values
    sum_hit_fake = abs(met_id_fc[hit_mask]).sum()
    score = np.empty(len(da_per))
    if sum_hit_fake == 0:  
        score[hit_mask] = np.nan
        return np.nan
    else:
        score[hit_mask] = abs(met_id_fc[hit_mask]) / sum_hit_fake
    score[~hit_mask] = -1 / (met_id_n_per - overlap_n_per)
     
    es_values = np.cumsum(score)
    max_index = np.abs(es_values).argmax()
    es_fake = es_values[max_index]
    return es_fake


def msea_ER(term_met, met_id, met_id_fc, forplot=False):
    term_met = term_met.split(';')
    overlap = list(set(term_met) & set(met_id))
    
    permutation_times = 1000
    
    if len(overlap) == 0:
        if forplot:
            return None
        else:
            return [0, None, None, None, None] + [None] * permutation_times
    
    da = pd.DataFrame({'met_id': met_id, 'met_id_fc': met_id_fc})
    da['met_id_fc'] = da['met_id_fc'].astype(float)
    da = da.sort_values(by='met_id_fc', ascending=False).reset_index(drop=True)
    da['rank'] = da.index + 1
    da['hit'] = np.where(da['met_id'].isin(overlap), 'yes', 'no')
    
    sum_hit = abs(da.loc[da['hit'] == 'yes', 'met_id_fc']).sum()
    if sum_hit == 0:
        da.loc[da['hit'] == 'yes', 'score'] = 0
        da.loc[da['hit'] == 'no', 'score'] = 0
        da['ES'] = da['score'].cumsum()
        es = 0
    else:
        da.loc[da['hit'] == 'yes', 'score'] = abs(da.loc[da['hit'] == 'yes', 'met_id_fc']) / sum_hit
        if(len(met_id) - len(overlap)) > 0:
            da.loc[da['hit'] == 'no', 'score'] = -1 / (len(met_id) - len(overlap))
        else:
            return [len(overlap), ';'.join(overlap), 0, 1, 0] + [np.nan]*permutation_times
        da['ES'] = da['score'].cumsum()
        es = da.loc[da['ES'].abs().idxmax(), 'ES']
        
    es_permu = []
    for i in range(permutation_times):
        es_permu.append(calculate_ES_fast((i, da, len(overlap), len(met_id))))
    es_permu = np.array(es_permu)

    es_permu = np.array(es_permu)
    pos_es_mean = es_permu[es_permu >= 0].mean() if len(es_permu[es_permu >= 0]) > 0 else 1
    neg_es_mean = es_permu[es_permu < 0].mean() if len(es_permu[es_permu < 0]) > 0 else 1
    
    nes_permu = np.where(es_permu >= 0, es_permu / pos_es_mean, es_permu / abs(neg_es_mean))
    
    # multiple ES values
    es = np.unique(es)
    if len(es) > 1:
        p_candidate = []
        nes_candidate = []
        for es_candidate in es:
            if es_candidate > 0:
                p_temp = np.sum(es_permu >= es_candidate) / np.sum(es_permu > 0)
                p_candidate.append(p_temp)
                nes_temp = es_candidate / pos_es_mean
                nes_candidate.append(nes_temp)
            elif es_candidate < 0:
                p_temp = np.sum(es_permu <= es_candidate) / np.sum(es_permu < 0)
                p_candidate.append(p_temp)
                nes_temp = es_candidate / abs(neg_es_mean)
                nes_candidate.append(nes_temp)
            else:
                p_candidate.append(1)
                nes_candidate.append(0)
        
        es_byp = es[np.argmin(p_candidate)]
        if len(es_byp) > 1:
            es_bynes = es[np.argmax(np.abs(nes_candidate))]
            if len(es_bynes) > 1:
                es = np.unique(np.abs(es))
            else:
                es = es_bynes
        else:
            es = es_byp
    
    if es >= 0:
        nes = es / pos_es_mean
        p_value = (es_permu >= es).sum() / (es_permu >= 0).sum()
    else:
        nes = es / abs(neg_es_mean)
        p_value = (es_permu <= es).sum() / (es_permu < 0).sum()
    
    if forplot:
        return da
    else:
        return [len(overlap), ';'.join(overlap), es, p_value, nes] + nes_permu.tolist()


def solve_redundant(hmdb_list, hmdb_list_value, choice="mean"):
    da = pd.DataFrame({"hmdb_list": hmdb_list, "hmdb_list_value": hmdb_list_value})
    
    if len(da["hmdb_list"].unique()) < len(da["hmdb_list"]):
        a = da["hmdb_list"].value_counts().reset_index()
        a.columns = ["hmdb_list", "Freq"]
        a = a[a["Freq"] > 1]
        
        if choice == "mean":
            for i in a["hmdb_list"]:
                da.loc[da["hmdb_list"] == i, "hmdb_list_value"] = da.loc[da["hmdb_list"] == i, "hmdb_list_value"].mean()
        elif choice == "max":
            for i in a["hmdb_list"]:
                da.loc[da["hmdb_list"] == i, "hmdb_list_value"] = da.loc[da["hmdb_list"] == i, "hmdb_list_value"].max()
        elif choice == "min":
            for i in a["hmdb_list"]:
                da.loc[da["hmdb_list"] == i, "hmdb_list_value"] = da.loc[da["hmdb_list"] == i, "hmdb_list_value"].min()
        da = da.drop_duplicates()
    return da


# main functon in MSEA
# the above functions are used in the main function

def msea_database(hmdb_list, hmdb_list_value, database, min_met=10, max_met=1000, choice="mean"):
    if database == "smpdb":
        res = smpdb_term_anno_new_MSEA
    elif database == "kegg":
        res = kegg_term_anno_new
    elif database == "humanGEM":
        res = humanGEM_term_new
    elif database == "reactome":
        res = reactome_term_anno_new
    elif database == "disease":
        res = disease_term_new
    elif database in ["Drug", "Disease", "Class", "Pathway"]:
        res = HUBMet_term[HUBMet_term['TermClass'] == database].copy()
        res.rename(columns={"HBM_ID": "HMDB.ID"}, inplace=True)
    
    if database in ["Drug", "Disease", "Class", "Pathway"]:
        fn_label = f"HUBMet_{database}"
    else:
        fn_label = database
            
    min_met_bp = res['N_met'].min()
    max_met_bp = res['N_met'].max()
    # Filter dataset
    res = res[(res['N_met'] >= min_met) & (res['N_met'] <= max_met)].copy()
    if len(res) == 0:
        print(f"No terms remained in {database} after applying the number filter.")
        print(f"The range of metabolite numbers in terms from {database} is from {min_met_bp} to {max_met_bp}.")
        return 0
        
    # redundant HMDB IDs
    temp = solve_redundant(hmdb_list, hmdb_list_value, choice)
    
    # parallel computing for each row (each row is a metabolite set)
    def process_row(hmdb_id, hmdb_list, hmdb_list_value):
        return msea_ER(hmdb_id, hmdb_list, hmdb_list_value)
    
    def parallel_msea(res, hmdb_list, hmdb_list_value):
        da_temp_need = Parallel(n_jobs=-1)(
            delayed(process_row)(row['HMDB.ID'], hmdb_list, hmdb_list_value)for _, row in res.iterrows())
        return da_temp_need
    da_temp_need = parallel_msea(res, temp['hmdb_list'], temp['hmdb_list_value'])
    
    da_temp_need = pd.DataFrame(da_temp_need)
    
    da_temp_need[2] = da_temp_need[2].apply(lambda x: x[0] if isinstance(x, np.ndarray) else x)
    da_temp_need[4] = da_temp_need[4].apply(lambda x: x[0] if isinstance(x, np.ndarray) else x)

    res.loc[:, ['hit_n', 'hit_HMDB_ID', 'ES', 'pvalue_nominal', 'NES']] = da_temp_need.iloc[:, :5].values
    if (res['hit_n'] == 0).all():
        print(f"No terms were mapped in {database}.")
        return 0
    NES_permutation_allset = da_temp_need.iloc[:, 5:]
    NES_permutation_allset = NES_permutation_allset.values.flatten()
    NES_permutation_allset = NES_permutation_allset[~np.isnan(NES_permutation_allset)]
    NES_permutation_allset = pd.to_numeric(pd.Series(NES_permutation_allset)).values
        
    # no overlap
    if res['ES'].isna().all():
        res.loc[:, 'FDR_pri'] = np.nan
        res.loc[:, 'FDR'] = np.nan
        return res[res['ES'].notna()]
    
    # Remove rows with no hits
    res = res[res['ES'].notna()]
    res['NES'] = res['NES'].astype(float)
    
    # Calculate FDR
    def calculate_fdr(x):
        if x >= 0:
            return (sum(NES_permutation_allset >= x) / sum(NES_permutation_allset >= 0)) / (sum(res['NES'] >= x) / sum(res['NES'] >= 0))
        else:
            return (sum(NES_permutation_allset <= x) / sum(NES_permutation_allset < 0)) / (sum(res['NES'] <= x) / sum(res['NES'] < 0))
    
    def parallel_calculate_fdr(values):
        results = Parallel(n_jobs=-1)(delayed(calculate_fdr)(x) for x in values)
        return results
    res['FDR_pri'] = parallel_calculate_fdr(res['NES'])

       
    # Sort by FDR
    res[['hit_n', 'N_met']] = res[['hit_n', 'N_met']].astype(int)
    col_na = ["ES", "pvalue_nominal", "NES", "FDR_pri"]
    for col in col_na:
        res[col] = res[col].astype(float)
    
    res = res.sort_values(by='FDR_pri', ascending=True)
    res['FDR'] = res['FDR_pri']
    res.loc[res['FDR_pri'] >= 1, 'FDR'] = 1
    
    table_export(res, database)
    res.to_csv(f"../Output/M2_MSEA/{fn_label}_MSEA_forRvisu.txt", sep='\t', index=False)
    return res




# test for HUBMet
testda_HBM = pd.read_csv('../testda/testda_hubmet_0709.txt', delimiter='\t')

 
import time
if __name__ == "__main__":
    for db in ["Class","Pathway","Disease","Drug"]:
        start_time = time.time()
        testda_msea_hbm = msea_database(hmdb_list=testda_HBM['metID'],
                                        hmdb_list_value=testda_HBM['fc'],
                                        database=db,
                                        min_met=10, max_met=1000)
        end_time = time.time()
        elapsed_time = end_time - start_time
        hours, rem = divmod(elapsed_time, 3600)
        minutes, seconds = divmod(rem, 60)
        formatted_time = "{:0>2}:{:0>2}:{:06.3f}".format(int(hours), int(minutes), seconds)
        print(f"Elapsed time: {formatted_time}")
