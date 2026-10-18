# How to publish this as a new Zenodo version

Do this AFTER the R2 manuscript is submitted in Editorial Manager.

1. Open the existing record: https://doi.org/10.5281/zenodo.19389356
2. Click **New version** (top right). This preserves the concept DOI and mints a new version DOI.
3. Remove the previous version's files if you want a clean version, or keep them; then upload the
   contents of this bundle: the whole `data/` and `scripts/` folders (drag the files in).
4. Set the metadata from `ZENODO_METADATA.txt` (title, version 1.3, license CC BY 4.0, keywords, description).
5. Under **Related/alternate identifiers**, confirm `isVersionOf 10.5281/zenodo.19389356` is present.
6. Save, then **Publish**. Publishing mints the version DOI and cannot be undone (you can only add further versions).
7. The manuscript cites the concept DOI, which does not change, so no manuscript edit is needed after publishing.

Note: publishing a Zenodo version is public and permanent. Review the file list and metadata before you click Publish.
