import numpy as np
import itertools as it
import pandas as pd

from typing import Union
from joblib import delayed, Parallel
from sklearn.model_selection import BaseCrossValidator, GridSearchCV
from sklearn.pipeline import FeatureUnion
from sklearn.base import clone as sk_clone

from rosvm.ranksvm.rank_svm_cls import KernelRankSVC, Labels
from rosvm.ranksvm.kernel_utils import generalized_tanimoto_kernel
from rosvm.feature_extraction.featurizer_cls import CircularFPFeaturizer, EStateIndFeaturizer


def run_baseline(rt_data: pd.DataFrame, estimator: Union[KernelRankSVC, GridSearchCV], cv: BaseCrossValidator,
                 featurizer: Union[FeatureUnion, CircularFPFeaturizer, EStateIndFeaturizer]):
    """
    """
    assert isinstance(estimator, GridSearchCV) or isinstance(estimator, KernelRankSVC)

    # Get list of all datasets
    dss = np.unique(rt_data["dataset"])

    # Run evaluation for each dataset as target
    scores = []
    for ds in dss:
        rt_data_ds = rt_data[rt_data["dataset"] == ds]

        # Extract features
        mols = rt_data_ds["smiles_iso"].values
        feat = sk_clone(featurizer).fit(mols)
        X = feat.transform(mols)

        # Set up label object
        y = Labels(rts=rt_data_ds["rt"].values, dss=ds)

        # Run the evaluation for each (train, test)-split
        for train, test in cv.split(X, y, groups=rt_data_ds["inchikey1"].values):
            est = sk_clone(estimator).fit(X[train], y[train])
            if isinstance(est, GridSearchCV):
                C_opt = est.best_estimator_.C
            else:
                C_opt = est.C

            scores.append([ds, est.score(X[test], y[test]), C_opt, len(feat)])

    return pd.DataFrame(scores, columns=["dataset", "cindex", "ranksvm_C", "n_fps"]).assign(setting="baseline")


def run_multisys(rt_data: pd.DataFrame, estimator: Union[KernelRankSVC, GridSearchCV], cv: BaseCrossValidator,
                 featurizer: CircularFPFeaturizer, local_fps_features: bool, use_classyfire_features: bool):
    """
    """
    assert isinstance(estimator, GridSearchCV) or isinstance(estimator, KernelRankSVC)

    # Get list of all datasets
    dss = np.unique(rt_data["dataset"])

    # Extract features
    if local_fps_features:
        n_fps = 0
        feat = FeatureUnion([(ds, sk_clone(featurizer)) for ds in dss])
        for ds, feat_ds in feat.transformer_list:
            n_fps += len(feat_ds.fit(rt_data[rt_data["dataset"] == ds]["smiles_iso"].values))
    else:
        feat = featurizer.fit(rt_data["smiles_iso"].values)
        n_fps = len(feat)
    K = generalized_tanimoto_kernel(feat.transform(rt_data["smiles_iso"].values))

    if isinstance(use_classyfire_features, str):
        assert use_classyfire_features in ["superclass", "class", "subclass"]
        Z, __dss, __cnts = get_classyfire_features(rt_data, use_classyfire_features, n_jobs=4)
        assert np.all(__dss == dss)
        K_Z = np.repeat(np.repeat(minmax_kernel_from_dict_lists(Z, n_jobs=4), __cnts, axis=0), __cnts, axis=1)
        K = K * K_Z

    # Set up label object
    y = Labels(rts=rt_data["rt"].values, dss=rt_data["dataset"].values)

    mols = rt_data["inchikey1"].values

    # Run evaluation for each dataset as target
    scores = []
    for ds in dss:
        idc_ds = y.get_idc_for_ds(ds)
        mols_ds = rt_data.iloc[idc_ds]["inchikey1"].values
        K_ds = K[idc_ds]  # shape: (n_samples_ds, n_samples_total)
        y_ds = y[idc_ds]

        # Run the evaluation for each (train, test)-split
        for _, test in cv.split(K_ds, y_ds, groups=mols_ds):
            train = [i for i in range(len(K)) if mols[i] not in mols_ds[test]]
            assert not any([m in mols[train] for m in mols_ds[test]])

            est = sk_clone(estimator).fit(K[np.ix_(train, train)], y[train])
            if isinstance(est, GridSearchCV):
                C_opt = est.best_estimator_.C
            else:
                C_opt = est.C

            scores.append([ds, est.score(K_ds[np.ix_(test, train)], y_ds[test]), C_opt])

    return pd.DataFrame(scores, columns=["dataset", "cindex", "ranksvm_C"]) \
        .assign(setting="multisys", n_fps=n_fps, use_classyfire_features=use_classyfire_features,
                local_fps_features=local_fps_features)


def run_leavesysout(rt_data: pd.DataFrame, estimator: Union[KernelRankSVC, GridSearchCV], cv: BaseCrossValidator,
                    featurizer: CircularFPFeaturizer, local_fps_features: bool, use_classyfire_features: bool):
    """
    """
    assert isinstance(estimator, GridSearchCV) or isinstance(estimator, KernelRankSVC)

    # Get list of all datasets
    dss = np.unique(rt_data["dataset"])

    # Extract features
    if local_fps_features:
        n_fps = 0
        feat = FeatureUnion([(ds, sk_clone(featurizer)) for ds in dss])
        for ds, feat_ds in feat.transformer_list:
            n_fps += len(feat_ds.fit(rt_data[rt_data["dataset"] == ds]["smiles_iso"].values))
    else:
        feat = featurizer.fit(rt_data["smiles_iso"].values)
        n_fps = len(feat)
    K = generalized_tanimoto_kernel(feat.transform(rt_data["smiles_iso"].values))

    if isinstance(use_classyfire_features, str):
        assert use_classyfire_features in ["superclass", "class", "subclass"]
        Z, __dss, __cnts = get_classyfire_features(rt_data, use_classyfire_features, n_jobs=4)
        assert np.all(__dss == dss)
        K_Z = np.repeat(np.repeat(minmax_kernel_from_dict_lists(Z, n_jobs=4), __cnts, axis=0), __cnts, axis=1)
        K = K * K_Z

    # Set up label object
    y = Labels(rts=rt_data["rt"].values, dss=rt_data["dataset"].values)

    mols = rt_data["inchikey1"].values

    # Run evaluation for each dataset as target
    scores = []
    for ds in dss:
        test = y.get_idc_for_ds(ds)
        mols_ds = rt_data.iloc[test]["inchikey1"].values
        K_ds = K[test]  # shape: (n_samples_ds, n_samples_total)
        y_ds = y[test]

        # Run the evaluation for each (train, test)-split
        train = [i for i in range(len(K)) if mols[i] not in mols_ds]
        assert not any([m in mols[train] for m in mols_ds])
        print("%s: n_train=%d, n_test=%d" % (ds, len(train), len(test)))

        est = sk_clone(estimator).fit(K[np.ix_(train, train)], y[train])
        if isinstance(est, GridSearchCV):
            C_opt = est.best_estimator_.C
        else:
            C_opt = est.C

        scores.append([ds, est.score(K_ds[:, train], y_ds), C_opt])

    return pd.DataFrame(scores, columns=["dataset", "cindex", "ranksvm_C"]) \
        .assign(setting="multisys", n_fps=n_fps, use_classyfire_features=use_classyfire_features,
                local_fps_features=local_fps_features)


def minmax_kernel_from_dicts(d1, d2):
    """
    A MinMax kernel implementation using two dictionaries with counts as input.
    """
    min_k = 0
    max_k = 0

    k_union = set(d1.keys()) | set(d2.keys())

    for key in k_union:
        v1 = d1.get(key, 0)
        v2 = d2.get(key, 0)

        min_k += np.minimum(v1, v2)
        max_k += np.maximum(v1, v2)

    return np.sum(min_k) / np.sum(max_k)


def minmax_kernel_from_dict_lists(l_d1, l_d2=None, n_jobs=4):
    if l_d2 is None:
        is_symmetric = True
        o_shape = (len(l_d1), len(l_d1))
    else:
        is_symmetric = False
        o_shape = (len(l_d1), len(l_d2))

    K = np.zeros(o_shape)

    if is_symmetric:
        res = Parallel(n_jobs=n_jobs)(delayed(minmax_kernel_from_dicts)(d1, d2)
                                      for d1, d2 in it.combinations(l_d1, 2))
        K[np.triu_indices(K.shape[0], k=1)] = res
        K[np.diag_indices(K.shape[0])] = 0.5
        K = K + K.T  # make symmetric

        assert np.all(np.equal(K[np.diag_indices(K.shape[0])], 1.0))
    else:
        res = Parallel(n_jobs=n_jobs)(delayed(minmax_kernel_from_dicts)(d1, d2)
                                      for d1, d2 in it.product(l_d1, l_d2))
        K = np.array(res).reshape(o_shape)

    return K


def _get_cls_cnts(l, normalize_counts):
    classes, counts = np.unique(l, return_counts=True)
    if normalize_counts:
        counts = np.atleast_1d(counts / np.sum(counts))
    return dict(zip(classes, counts))


def get_classyfire_features(rt_data, level, n_jobs=4, normalize_counts=True):
    dss, cts = np.unique(rt_data.dataset, return_counts=True)
    feat = Parallel(n_jobs=n_jobs)(delayed(_get_cls_cnts)(rt_data[rt_data.dataset == ds][level], normalize_counts)
                                   for ds in dss)
    return feat, dss, cts


if __name__ == "__main__":
    from sklearn.model_selection import GridSearchCV, KFold, GroupShuffleSplit
    from rosvm.ranksvm.rank_svm_cls import KernelRankSVC
    from rosvm.feature_extraction.featurizer_cls import CircularFPFeaturizer

    rt_data = pd.read_csv("../metfrag/mol_rt_info/all_rt_data.csv").sort_values(by="dataset")

    run_leavesysout(
        rt_data,
        cv=GroupShuffleSplit(test_size=0.25, n_splits=25, random_state=33),
        estimator=GridSearchCV(KernelRankSVC(kernel="precomputed", random_state=932),
                               param_grid={"C": [1/8, 0.25, 0.5, 1, 4, 8, 16, 32]}, cv=KFold(n_splits=10), n_jobs=4),
        featurizer=CircularFPFeaturizer(only_freq_subs=True, min_subs_freq=0.05, output_dense_matrix=True, radius=6,
                                        use_chirality=True),
        use_classyfire_features=False, local_fps_features=False)
