
import os
import numpy as np

YANK_LIGANDS = {}

YANK_LIGANDS["1-methylpyrrole.A__AAA"]      = "methylpyrrole"
YANK_LIGANDS["benzene.A__AAA"]              = "benzene"
YANK_LIGANDS["lysozyme.active.A__AAK"]      = "benzofuran"
YANK_LIGANDS["lysozyme.active.A__AAZ"]      = "allyl ethyl sulfide"
YANK_LIGANDS["lysozyme.active.A__ABA"]      = "hexafluorobenzene"
YANK_LIGANDS["lysozyme.active.A__ABC"]      = "indole"
YANK_LIGANDS["lysozyme.active.A__ABG"]      = "m-xylene"
YANK_LIGANDS["lysozyme.active.A__ABJ"]      = "n-hexylbenzene"
YANK_LIGANDS["lysozyme.active.A__ABK"]      = "nitrobenzene"
YANK_LIGANDS["lysozyme.active.A__AAU"]      = "4-ethyltoluene"

YANK_LIGANDS["lysozyme.inactive.A__AA0"]    = "ethanol"
YANK_LIGANDS["lysozyme.inactive.A__AA1"]    = "methanol"
YANK_LIGANDS["lysozyme.inactive.A__AA3"]    = "dimethyl-sulfoxide"
YANK_LIGANDS["lysozyme.inactive.A__AAS"]    = "DL-camphor"
YANK_LIGANDS["lysozyme.inactive.A__AAY"]    = "1-propanol"
YANK_LIGANDS["lysozyme.inactive.A__AAZ"]    = "1,1-diethylurea"
YANK_LIGANDS["lysozyme.inactive.A__ABA"]    = "1,4 diiodobenzene"
YANK_LIGANDS["lysozyme.inactive.A__ABB"]    = "1,2-diiodobenzene"
YANK_LIGANDS["lysozyme.inactive.A__ABH"]    = "2-bromoethanol"
YANK_LIGANDS["lysozyme.inactive.A__ABI"]    = "2-iodoethanol"
YANK_LIGANDS["lysozyme.inactive.A__ABJ"]    = "benzyl alcohol"
YANK_LIGANDS["lysozyme.inactive.A__ABK"]    = "benzaldehyde oxime"
YANK_LIGANDS["phenol.A__AAA"]               = "phenol"

YANK_LIGANDS["p-xylene.A__AAA"]             = "p-xylene"


SIX_REF_SYSTEMS = ["lysozyme.active.A__ABJ", "benzene.A__AAA", "lysozyme.inactive.A__AAS",
                    "1-methylpyrrole.A__AAA", "phenol.A__AAA", "p-xylene.A__AAA"]


def load_scores(file, id_col, score_col, std_col, exclude_ligands):
    scores = {}
    standard_devs = {}

    with open(file, "r") as handle:
        for line in handle:
            if not line.strip().startswith("#"):
                entries = line.split()
                id = entries[id_col]
                val = entries[score_col]
                std = entries[std_col]

                if id not in exclude_ligands:
                    if val.lower() not in ["inf", "nan"]:
                        scores[id] = np.float(val)
                        standard_devs[id] = np.float(std)
    return scores, standard_devs


def matching_scores(score1, score2, err1, err2, allowed_ligands=[]):
    ligands = set(score1.keys()).intersection(set(score2.keys()))

    if len(allowed_ligands) > 0:
        ligands = ligands.intersection(allowed_ligands)
    ligands = list(ligands)

    x, y, xerr, yerr = [], [], [], []

    for ligand in ligands:
        x.append(score1[ligand])
        y.append(score2[ligand])

        xerr.append(err1[ligand])
        yerr.append(err2[ligand])

    return np.array(x), np.array(y), np.array(xerr), np.array(yerr), ligands


def write_pairs(score1, score2, err1, err2, out, allowed_ligands):
    ligands = set(score1.keys()).intersection( set(score2.keys()))
    if len(allowed_ligands) > 0:
        ligands = ligands.intersection( set(allowed_ligands))
    ligands = list(ligands)

    diff = {l:np.abs(score1[l] - score2[l]) for l in ligands}
    ligands.sort(key=lambda l: diff[l], reverse=True)

    with open(out, "w") as handle:
        handle.write("# ligand     x     y    xerr    yerr   diff\n")
        for l in ligands:
            handle.write("%30s %20.10f %20.10f %20.10f %20.10f %20.10f\n" % (l, score1[l], score2[l],
                                                                             err1[l], err2[l], diff[l]))
    return None


def snapshot_indices_for_algdock(nframes_per_state=10000):
    """
    :param nframes_per_state: int
    :return: [5000, 9000]
             index of snaphots from each state chosen for algdock calculations
    """
    first_stride = 1000
    first_list = [i for i in range(nframes_per_state) if i%first_stride == 0]

    last_list = [first_list[5+10*i] for i in range(len(first_list)) if 5+10*i < len(first_list) ]
    last_list.extend( [ first_list[9+10*i] for i in range(len(first_list)) if 9+10*i < len(first_list) ] )
    return sorted(last_list)


def load_interaction_energies(snaphot_col=0, value_col=1, path="OpenMM_OBC2_interaction_energies"):
    """
    :param snaphot_col: int
    :param value_col: int
    :param path: str
    :return: interaction_energies, dict
            interaction_energies[ligand][snapshot] -> float
    """
    interaction_energies = {}
    for ligand in SIX_REF_SYSTEMS:
        filename = os.path.join(path, ligand + ".dat")
        print("loaing " + filename)

        interaction_energies[ligand] = {}
        with open(filename, "r") as handle:
            for line in handle:
                if "#" not in line:
                    entries = line.split()
                    snapshot = entries[snaphot_col]
                    value = entries[value_col]
                    interaction_energies[ligand][snapshot] = np.float(value)
    return interaction_energies


def load_bpmfs(file_name, exclude_nan=False):
    """
    :param file_name: str
    :param exclude_nan: bool, whether to exclude NaN values
    :return: bpmfs, dict
             bpmfs[snapshot]  -> float
             e.g. bpmfs["4"]  -> -2.8192
    """
    bpmfs = {}
    with open(file_name, "r") as handle:
        for line in handle:
            snapshot, value = line.split()

            if value.lower() != "nan":
                bpmfs[snapshot] = np.float(value)
            else:
                if not exclude_nan:
                    bpmfs[snapshot] = np.nan
    return bpmfs