# smilite

smilite is a Python module to download and analyze SMILES strings (Simplified Molecular-Input Line-entry System) of chemical compounds from ZINC (a free database of commercially-available compounds for virtual screening, [http://zinc.docking.org](http://zinc.docking.org)).  
Now supports both Python 3.x and Python 2.x.


![](https://raw.github.com/rasbt/smilite/master/images/smilite_overview.png)  

#### Sections

• [Installation](#installation)  
• [Simple command line online query scripts](#simple_cmd_scripts)  
      - [lookup_zincid.py](#lookup_zincid)  
      - [lookup_smile_str.py](#lookup_smile_str)  
• [CSV file command line scripts](#csv_scripts)  
      - [gen_zincid_smile_csv.py (downloading SMILES)](#gen_zincid)  
      - [comp_smile_strings.py (checking for duplicates within 1 file)](#comp_smile)  
      - [comp_2_smile_files.py (checking for duplicates across 2 files)](#comp_2_smile)  
• [SQLite file command line scripts](#sqlite_scripts)  
      - [lookup_single_id.py](#lookup1id)  
      - [lookup_smile.py](#lookupsmile)  
      - [add_to_sqlite.py](#add_to_sqlite)  
      - [sqlite_to_csv.py](#sqlite_to_csv)  
• [Changelog](#changelog)  

<a name="installation"></a>

# Installation

You can use the following command to install smilite:  
`pip install smilite`  
or  
`easy_install smilite`

Alternatively, you can download the package manually from the Python Package Index [https://pypi.python.org/pypi/smilite](https://pypi.python.org/pypi/smilite), unzip it, navigate into the package, and use the command:

`python3 setup.py install`

<a name="simple_cmd_scripts"></a>

# Simple command line online query scripts

If you downloaded the smilite package from [https://pypi.python.org/pypi/smilite](https://pypi.python.org/pypi/smilite) or [https://github.com/rasbt/smilite](https://github.com/rasbt/smilite), you can use the command line scripts I provide in the `scripts/cmd_line_online_query_scripts` dir.

<a name="lookup_zincid"></a>

### lookup_zincid.py

Retrieves the SMILES string and simplified SMILES string for a given ZINC ID  
from the online Zinc. It uses [ZINC12](http://zinc.docking.org) as the default backend, and via an additional commandline argument `zinc15`, the [ZINC15](http://zinc15.docking.org) database will be used instead.

**Usage:**  
`[shell]>> python3 lookup_zincid.py ZINC_ID [zinc12/zinc15]`  

**Example (retrieve data from ZINC):**  
`[shell]>> python3 lookup_zincid.py ZINC01234567 zinc15`  

**Output example:**

<pre>ZINC01234567
C[C@H]1CCCC[NH+]1CC#CC(c2ccccc2)(c3ccccc3)O
CC1CCCCN1CCCC(C2CCCCC2)(C3CCCCC3)O
</pre>

Where  
- 1st row: ZINC ID  
- 2nd row: SMILES string  
- 3rd row: simplified SMILES string

<a name="lookup_smile_str"></a>

### lookup_smile_str.py

Retrieves the corresponding ZINC_IDs for a given SMILES string  
from the online ZINC database. 

**Usage:**  
`[shell]>> python3 lookup_smile_str.py SMILE_str`  

**Example (retrieve data from ZINC):**  
`[shell]>> python3 lookup_smile_str.py "C[C@H]1CCCC[NH+]1CC#CC(c2ccccc2)(c3ccccc3)O"`  

**Output example:**

<pre>ZINC01234567
ZINC01234568
ZINC01242053
ZINC01242055</pre>

<a name="csv_scripts"></a>

# CSV file command line scripts

If you downloaded the smilite package from [https://pypi.python.org/pypi/smilite](https://pypi.python.org/pypi/smilite) or [https://github.com/rasbt/smilite](https://github.com/rasbt/smilite), you can use the command line scripts I provide in the `scripts/csv_scripts` dir.

<a name="gen_zincid"></a>

### gen_zincid_smile_csv.py (downloading SMILES)

Generates a ZINC_ID,SMILE_STR csv file from a input file of ZINC IDs. The input file should consist of 1 columns with 1 ZINC ID per row. [ZINC12](http://zinc.docking.org) is used as the default backend, and via an additional commandline argument `zinc15`, the [ZINC15](http://zinc15.docking.org) database can be used instead.

**Usage:**  
`[shell]>> python3 gen_zincid_smile_csv.py in.csv out.csv [zinc12/zinc15]`

**Example:**  
`[shell]>> python3 gen_zincid_smile_csv.py ../examples/zinc_ids.csv ../examples/zid_smiles.csv zinc15`

**Screen Output:**

<pre>Downloading SMILES
0%                          100%
[##########                    ] | ETA[sec]: 106.525 </pre>

**Input example file format:**  
![](https://raw.github.com/rasbt/smilite/master/images/zinc_ids.png)  
[zinc_ids.csv](https://raw.github.com/rasbt/smilite/master/examples/zinc_ids.csv)

**Output example file format:**  
![](https://raw.github.com/rasbt/smilite/master/images/zid_smiles.png)  
[zid_smiles.csv](https://raw.github.com/rasbt/smilite/master/examples/zid_smiles.csv)

<a name="comp_smile"></a>

### comp_smile_strings.py (checking for duplicates within 1 file)

Compares SMILES strings within a 2 column CSV file (ZINC_ID,SMILE_string) to identify duplicates. Generates a new CSV file with ZINC IDs of identified duplicates listed in a 3rd-nth column(s).

**Usage:**  
`[shell]>> python3 comp_smile_strings.py in.csv out.csv [simplify]`

**Example 1:**  
`[shell]>> python3 comp_smile_strings.py ../examples/zinc_smiles.csv ../examples/comp_smiles.csv`

**Input example file format:**  
![](https://raw.github.com/rasbt/smilite/master/images/zid_smiles.png)  
[zid_smiles.csv](https://raw.github.com/rasbt/smilite/master/examples/zid_smiles.csv)

**Output example file format 1:**  
![](https://raw.github.com/rasbt/smilite/master/images/comp_smiles.png)  
[comp_smiles.csv](https://raw.github.com/rasbt/smilite/master/examples/comp_smiles.csv)

Where  
- 1st column: ZINC ID  
- 2nd column: SMILES string  
- 3rd column: number of duplicates  
- 4th-nth column: ZINC IDs of duplicates

**Example 2:**  
`[shell]>> python3 comp_smile_strings.py ../examples/zid_smiles.csv ../examples/comp_simple_smiles.csv simplify`

**Output example file format 2:** ![](https://raw.github.com/rasbt/smilite/master/images/comp_simple_smiles.png)  
[comp_simple_smiles.csv](https://raw.github.com/rasbt/smilite/master/examples/comp_simple_smiles.csv)

<a name="comp_2_smile"></a>

### comp_2_smile_files.py (checking for duplicates across 2 files)

Compares SMILES strings between 2 input CSV files, where each file consists of rows with 2 columns ZINC_ID,SMILE_string to identify duplicate SMILES string across both files.  
Generates a new CSV file with ZINC IDs of identified duplicates listed in a 3rd-nth column(s).

**Usage:**  
`[shell]>> python3 comp_2_smile_files.py in1.csv in2.csv out.csv [simplify]`

**Example:**  
`[shell]>> python3 comp_2_smile_files.py ../examples/zid_smiles2.csv ../examples/zid_smiles3.csv ../examples/comp_2_files.csv`

**Input example file 1:**  
![](https://raw.github.com/rasbt/smilite/master/images/zid_smiles2.png)  
[zid_smiles2.csv](https://raw.github.com/rasbt/smilite/master/examples/zid_smiles2.csv)

**Input example file 2:**  
![](https://raw.github.com/rasbt/smilite/master/images/zid_smiles3.png)  
[zid_smiles3.csv](https://raw.github.com/rasbt/smilite/master/examples/zid_smiles3.csv)

**Output example file format:**  
![](https://raw.github.com/rasbt/smilite/master/images/comp_2_files.png)  
[comp_2_files.csv](https://raw.github.com/rasbt/smilite/master/examples/comp_2_files.csv)

Where:  
- 1st column: name of the origin file  
- 2nd column: ZINC ID  
- 3rd column: SMILES string  
- 4th-nth column: ZINC IDs of duplicates

<a name="sqlite_scripts"></a>

# SQLite file command line scripts

If you downloaded the smilite package from [https://pypi.python.org/pypi/smilite](https://pypi.python.org/pypi/smilite) or [https://github.com/rasbt/smilite](https://github.com/rasbt/smilite), you can use the command line scripts I provide in the `scripts/sqlite_scripts` dir.

<a name="lookup1id"></a>

### lookup_single_id.py

Retrieves the SMILES string and simplified SMILES string for a given ZINC ID  
from a previously built smilite SQLite database or from the online ZINC database.

**Usage:**  
`[shell]>> python3 lookup_single_id.py ZINC_ID [sqlite_file]`  

**Example1 (retrieve data from a smilite SQLite database):**  
`[shell]>> python3 lookup_single_id.py ZINC01234567 ~/Desktop/smilite_db.sqlite`  

**Example2 (retrieve data from the ZINC online database):**  
`[shell]>> python3 lookup_single_id.py ZINC01234567`  

**Output example:**

<pre>ZINC01234567
C[C@H]1CCCC[NH+]1CC#CC(c2ccccc2)(c3ccccc3)O
CC1CCCCN1CCCC(C2CCCCC2)(C3CCCCC3)O
</pre>

Where  
- 1st row: ZINC ID  
- 2nd row: SMILES string  
- 3rd row: simplified SMILES string

<a name="lookupsmile"></a>

### lookup_smile.py

Retrieves the ZINC ID(s) for a given SMILES string or simplified SMILES string from a previously built smilite SQLite database.

**Usage:**  
`[shell]>> python3 lookup_smile.py sqlite_file SMILE_STRING [simplify]`  

**Example1 (search for SMILES string):**  
`[shell]>> python3 lookup_smile.py ~/Desktop/smilite.sqlite "C[C@H]1CCCC[NH+]1CC#CC(c2ccccc2)(c3ccccc3)O"`  

**Example2 (search for simplified SMILES string):**  
`[shell]>> python3 lookup_smile.py ~/Desktop/smilite.sqlite "CC1CCCCN1CCCC(C2CCCCC2)(C3CCCCC3)O" simple`  

**Output example:**

<pre>ZINC01234567
C[C@H]1CCCC[NH+]1CC#CC(c2ccccc2)(c3ccccc3)O
CC1CCCCN1CCCC(C2CCCCC2)(C3CCCCC3)O
</pre>

Where  
- 1st row: ZINC ID  
- 2nd row: SMILES string  
- 3rd row: simplified SMILES string

<a name="add_to_sqlite"></a>

### add_to_sqlite.py

Reads ZINC IDs from a CSV file and looks up SMILES strings and simplified SMILES strings from the ZINC online database. Writes those SMILES strings to a smilite SQLite database. A new database will be created if it doesn't exist, yet.

**Usage:**  
`[shell]>> python3 add_to_sqlite.py sqlite_file csv_file`  

**Example:**  
`[shell]>> python3 add_to_sqlite.py ~/Desktop/smilite.sqlite ~/Desktop/zinc_ids.csv`  

**Input CSV file example format:**

<pre>ZINC01234567
ZINC01234568
...
</pre>

An example of the smilite SQLite database contents after successful insertion is shown in the image below. ![https://raw.github.com/rasbt/smilite/master/images/add_to_sqlite_1.png](https://raw.github.com/rasbt/smilite/master/images/add_to_sqlite_1.png)

<a name="sqlite_to_csv"></a>

### sqlite_to_csv.py

Writes contents of an SQLite smilite database to a CSV file.

**Usage:**  
`[shell]>> python3 sqlite_to_csv.py sqlite_file csv_file`  

**Example:**  
`[shell]>> python3 sqlite_to_csv.py ~/Desktop/smilite.sqlite ~/Desktop/zinc_smiles.csv`  

**Input CSV file example format:**

<pre>ZINC_ID,SMILE,SIMPLE_SMILE
ZINC01234568,C[C@@H]1CCCC[NH+]1CC#CC(c2ccccc2)(c3ccccc3)O,CC1CCCCN1CCCC(C2CCCCC2)(C3CCCCC3)O
ZINC01234567,C[C@H]1CCCC[NH+]1CC#CC(c2ccccc2)(c3ccccc3)O,CC1CCCCN1CCCC(C2CCCCC2)(C3CCCCC3)O
</pre>

An example of the CSV file contents opened in an spreadsheet program is shown in the image below. ![https://raw.github.com/rasbt/smilite/master/images/sqlite_to_csv_2.png](https://raw.github.com/rasbt/smilite/master/images/sqlite_to_csv_2.png)



<a name="changelog"></a>

# Changelog

**VERSION 2.3.1 (07/25/2020)**

- Fix bug to allow `zinc15` option in gen_zincid_smile_csv.py script 

**VERSION 2.3.0 (06/10/2020)**

- Fixes ZINC URL in `lookup_smile_str.py`
- Adds an optional command line parameter (with arguments `zinc15` or `zinc12`) for `lookup_smile_str.py`

**VERSION 2.2.0**

*   Provides an optional command line argument (zinc15) to use ZINC15 as a backend for downloading SMILES

**VERSION 2.1.0**

*   Functions and scripts to fetch ZINC IDs corresponding to a SMILES string query

**VERSION 2.0.1**

*   Progress bar for add_to_sqlite.py

**VERSION 2.0.0**

*   added SQLite features

**VERSION 1.3.0**

*   added script and module function to compare SMILES strings across 2 files.

**VERSION 1.2.0**

*   added Python 2.x support

**VERSION 1.1.1**

*   PyPrind dependency fix

**VERSION 1.1.0**

*   added a progress bar (PyPrind) to `generate_zincid_smile_csv()` function