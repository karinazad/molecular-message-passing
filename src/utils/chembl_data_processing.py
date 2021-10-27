import os
import ast
from typing import Optional
import pandas as pd

from chembl_webresource_client.query_set import QuerySet


def get_bioactivity_compound_data(bioacitivity_query, compound_query, data_path, target_name):
    """
    ----------
    bioacitivity_query
    compound_query
    data_path
    target_name

    Returns
    -------

    """

    # Get and process bioactivity dataset from the provided query
    print("1. Get and process bioactivity dataset\n")
    bioactivities_df = retrieve_or_create_dataframe_from_query(
        query=bioacitivity_query,
        query_name="bioactivities",
        data_path=data_path,
        target_name=target_name,
    )

    # Remove NaNs and duplicates, convert measurement column to float, and keep only entries in standard units
    bioactivities_df = preprocess_queried_dataset(
        df=bioactivities_df,
        columns_to_float=["standard_value"],
        standard_unit="nM"
    )

    # Get and process compounds dataset from the provided query
    print("2. Get and process compounds dataset\n")
    compounds_df = retrieve_or_create_dataframe_from_query(
        query=compound_query,
        query_name="compounds",
        data_path=data_path,
        target_name=target_name,
    )

    # Remove NaNs and duplicates
    compounds_df = preprocess_queried_dataset(
        df=compounds_df,
    )

    # Extract SMILES representation from molecular representation of a compound
    compounds_df = extract_smiles_from_molecular_representation(
        df=compounds_df
    )

    # Merge bioactivity and compound dataset on ChemBL ID
    print("3. Merge bioactivity and compound dataset\n")
    merged_df = pd.merge(
        bioactivities_df[["molecule_chembl_id", "standard_value", "standard_units"]],
        compounds_df,
        on="molecule_chembl_id",
    )

    # Reset row indices
    merged_df.reset_index(drop=True, inplace=True)

    print(f"\tFinal number of samples:  {merged_df.shape[0]}")

    return merged_df


def preprocess_queried_dataset(
        df: pd.DataFrame,
        columns_to_float: Optional[list] = None,
        drop_na: bool = True,
        drop_duplicates: bool = True,
        reset_index=True,
        standard_unit: Optional[str] = None,
):
    """
    1. Convert standard_value's datatype from object to float
    2. Delete entries with missing values
    3. Keep only entries with standard_unit == nM
    4. Delete duplicate molecules
    5. Reset DataFrame index

    Parameters
    ----------
    df
    columns_to_float
    drop_na
    drop_duplicates
    reset_index
    standard_unit

    Returns
    -------

    """
    n_samples = len(df)
    print(f"\tInitial number of samples: {n_samples} \n")

    # Convert standard_value's datatype from object to float
    if columns_to_float is not None:
        df = df.astype({column: "float64" for column in columns_to_float})

    # Delete entries with missing values
    if drop_na:
        df.dropna(axis=0, how="any", inplace=True)
        print(f"\t\t Entries with NaNs: {n_samples - len(df)}."
              f"\n\t\t New number of samples: {len(df)} \n")
        n_samples = len(df)

    # 3 Keep only entries with standard_unit
    if standard_unit != "all" and standard_unit is not None:
        assert "standard_units" in df.columns

        df = df[df["standard_units"] == standard_unit]
        print(f"\t\t Entries with non-standard units: {n_samples - len(df)}"
              f"\n\t\t New number of samples: {len(df)}\n")
        n_samples = len(df)

    # Delete duplicate molecules
    if drop_duplicates:
        df.drop_duplicates("molecule_chembl_id", keep="first", inplace=True)
        print(f"\t\t Entries duplicates: {n_samples - len(df)}"
              f"\n\t\t New number of samples: {len(df)}\n")
        n_samples = len(df)

    # Reset DataFrame index
    if reset_index:
        df.reset_index(drop=True, inplace=True)

    print(f"\tNumber of samples: {n_samples} \n")

    return df


def retrieve_or_create_dataframe_from_query(
        query: QuerySet,
        data_path: str,
        query_name: str,
        target_name: str = "",
        save_df: bool = True,
):
    data_file_path = os.path.join(data_path, f"{target_name}_{query_name}.csv")

    # Try loading saved data
    try:
        df = pd.read_csv(data_file_path)
        df.drop(["Unnamed: 0"], inplace=True, axis=1)

    # Or load into a dataframe and save the datta
    except FileNotFoundError:
        df = pd.DataFrame.from_records(query)

        if save_df:
            if not os.path.exists(data_path):
                os.mkdir(data_path)

            df.to_csv(data_file_path)

    return df


def extract_smiles_from_molecular_representation(df):
    """

    Parameters
    ----------
    df

    Returns
    -------

    """
    assert "molecule_structures" in df.columns

    # Convert strings to dictionaries if needed
    if type(df.iloc[0].molecule_structures) is not dict:
        df.molecule_structures = df.molecule_structures.apply(lambda x: ast.literal_eval(x))

    # Keep only SMILES representation
    df["canonical_smiles"] = [
        compounds["molecule_structures"]["canonical_smiles"] for _, compounds in df.iterrows()
    ]

    df.drop("molecule_structures", axis=1, inplace=True)
    df.dropna(axis=0, how="any", inplace=True)

    return df
