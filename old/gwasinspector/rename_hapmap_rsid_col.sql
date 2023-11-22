-- There's a bug/typo in GWASinspector where it is expecting the rsID column in the SQLite file to be "REF_RSID", not "ID". σ_σ
-- This script is simply to rename the "ID" column to "REF_RSID" in HapMap_CEU_r28_b36_EDIT_v10c_v5.sqlite
-- Use this script like so: sqlite3 <filename> < this_script.sql

CREATE TABLE tmp_table (
  hID TEXT,
  REF_RSID TEXT,
  REF TEXT,
  ALT TEXT,
  AF REAL,
  TSA TEXT,
  ignore numeric
);

-- Copy data from existing "variants" table to tmp_table
INSERT INTO tmp_table (hID, REF_RSID, REF, ALT, AF, TSA, ignore)
  SELECT hID, ID AS REF_RSID, REF, ALT, AF, TSA, ignore
  FROM variants;

-- Drop the old table
DROP TABLE variants;

-- Rename new table to match the old name
ALTER TABLE tmp_table RENAME TO variants;
