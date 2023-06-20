import sys
import glob


def get_outdock_distribution(outdock_fn):
    with open(outdock_fn) as f:
        lines = f.readlines()
    i = 0

    data_dict = dict()
    found_open = False
    while i < len(lines):
        line = lines[i]
        if line.startswith(' open'):
            found_open = True
            break
        i += 1
    if not found_open:
        i = 0
    while i < len(lines):
        line = lines[i]
        if not (line.startswith(' open') or line.startswith(' close') or line.startswith(' /')):
            ll = line.strip().split()
            if len(ll) >= 2:
                if ll[1].startswith('NC') or ll[1].startswith('kezinc'): #'ZINC' is truncated to 'NC' in dock3.8
                    zincid = ll[1]
                    try:
                        score = float(ll[-1])
                    except ValueError:
                        i += 1
                        continue
                    if zincid in data_dict:
                        if score < data_dict[zincid]:
                            data_dict[zincid] = score
                    else:
                        data_dict[zincid] = score
        i += 1
    return data_dict


def fuse_data_dicts(ref_dict, add_dict):
    for key in add_dict:
        if key not in ref_dict:
            ref_dict[key] = add_dict[key]
        else:
            newval = add_dict[key]
            if newval < ref_dict[key]:
                ref_dict[key] = newval
    return ref_dict


#calls the other two methods
def write_scores_df(outname, rgx):
    #if rgx is None:
    #    rgx = f"/nfs/home/jborowsky/rotation/ampc-v1/run4/dock-output{index}/*/OUTDOCK.0" #'/nfs/home/jborowsky/rotation/d4-v1/run18/controls/controls*/OUTDOCK'
    #    #rgx = '/wynton/group/bks/work/omailhot/complement_receptors/c5ar_goldilocks/output/*/OUTDOCK.0'
    outdocks = glob.glob(rgx)
    dd = dict()
    for od in outdocks:
        d = get_outdock_distribution(od)
        #print(len(d), len(dd))
        dd = fuse_data_dicts(dd, d)
        #print('{} done'.format(od))
    
    ddsorted = sorted(dd.items(), key=lambda x:x[1])
        
    with open(outname, 'w') as f:
        f.write('zid score\n')
        for zid in dd:
            f.write('{} {}\n'.format(zid, dd[zid]))

    return ddsorted


def get_best_smiles(outname, smiles_path, rgx):
    scores_sorted = write_scores_df(outname, rgx)
    #print(scores_sorted)

    #print(f"Best molecule: {scores_sorted[0]}")

    best_trunc_zinc = scores_sorted[0][0].strip()[:-2]   

    fs = open(smiles_path, "r")
    for line in fs:
        
        if line.split(" ")[1].strip()[2:] == best_trunc_zinc:
            
            return [line.split(" ")[0], line.split(" ")[1], scores_sorted[0][1]]

#best_smiles = get_best_smiles("dock38-ampc-36k-analogues-scores", "/nfs/home/jborowsky/rotation/ampc-v1/input-files/input-smiles/ampc-exp-analogues-i1521-o36331.smi")

#print(best_smiles)
