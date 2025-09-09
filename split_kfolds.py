import os, sys
import anndata as ad

sys.path.insert(1, '/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Inkar/scripts/bmark')
import bmark
from bmark.utils.dataset import get_kfold_ind


# saves to current location + \folds + \fold_# + h5ads
import os
import anndata as ad
from tqdm.auto import tqdm   # instead of tqdm.auto

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)


def split_dataset(taxonomy, stratify_by, save_folds_dir, n_folds=10, show_progress=True):
    kfold_save_dir = os.path.join(save_folds_dir, "folds")
    os.makedirs(kfold_save_dir, exist_ok=True)

    taxonomy.uns["benchmark"] = {"k_fold": {}}

    pbar = tqdm(range(n_folds), disable=not show_progress, desc="Creating folds")
    for fold in pbar:
        train_ind, validation_ind = get_kfold_ind(
            obs=taxonomy.obs,
            stratify_by=stratify_by,
            fold=fold,
            n_folds=n_folds
        )
        taxonomy.uns["benchmark"]["k_fold"][f"fold_{fold+1}"] = {
            "train_ind": train_ind,
            "val_ind": validation_ind
        }

        train_taxonomy = taxonomy[train_ind, :]
        val_taxonomy   = taxonomy[validation_ind, :]

        fold_folder_path = os.path.join(kfold_save_dir, f'fold_{fold+1}')
        os.makedirs(fold_folder_path, exist_ok=True)

        # Optional: small inner bar for the two writes
        if show_progress:
            with tqdm(total=2, desc=f"Fold {fold+1} writes", leave=False) as wbar:
                train_taxonomy.write_h5ad(os.path.join(fold_folder_path, 'train.h5ad'))
                wbar.update(1)
                val_taxonomy.write_h5ad(os.path.join(fold_folder_path, 'validation.h5ad'))
                wbar.update(1)
        else:
            train_taxonomy.write_h5ad(os.path.join(fold_folder_path, 'train.h5ad'))
            val_taxonomy.write_h5ad(os.path.join(fold_folder_path, 'validation.h5ad'))

        # Nice live info in the main bar
        pbar.set_postfix(fold=fold+1, train=len(train_ind), val=len(validation_ind))
