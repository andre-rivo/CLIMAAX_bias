{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "01d9f030-ea3f-42fc-90fc-581954b51bdd",
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "776ba803-adf2-41a2-88ea-de6026bceb5e",
   "metadata": {
    "tags": []
   },
   "outputs": [],
   "source": [
    "import numpy as np\n",
    "import geopandas as gpd\n",
    "import os\n",
    "import xarray as xr\n",
    "import matplotlib.pyplot as plt\n",
    "import glob\n",
    "import rioxarray as rio\n",
    "from tqdm import tqdm\n",
    "import plotly.graph_objects as go\n",
    "import pandas as pd\n",
    "import plotly.express as px\n",
    "import plotly.colors as pcolors\n",
    "from dash import dcc, html, Input, Output, Dash"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "id": "3a6ed878-c693-4fe9-9b79-cc8c0837098e",
   "metadata": {
    "tags": []
   },
   "outputs": [],
   "source": [
    "#import eobs tas\n",
    "\n",
    "coords = ['longitude','latitude']\n",
    "eobs_tas = xr.open_dataset('tg_ens_mean_0.1deg_reg_v30.0e_REMAP_FAIZAN_mean1971-2000.nc')\n",
    "eobs_tas=eobs_tas.assign_coords({coord: eobs_tas[coord].round(3) for coord in coords})\n",
    "    "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "id": "30145dec-fdc1-413e-870c-08a248f3fdbf",
   "metadata": {
    "tags": []
   },
   "outputs": [],
   "source": [
    "#import eobs pr\n",
    "coords = ['longitude','latitude']\n",
    "eobs_pr = xr.open_dataset('rr_ens_mean_0.1deg_reg_v30.0e_REMAP_FAIZAN_MEAN_YEARSUM_1971-2000.nc')\n",
    "eobs_pr=eobs_pr.assign_coords({coord: eobs_pr[coord].round(3) for coord in coords})\n",
    "    "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "id": "634a7890-da8b-4848-8f2c-ce1dfb16d37c",
   "metadata": {
    "tags": []
   },
   "outputs": [],
   "source": [
    "gcm_list = ['EC-EARTH','HadGEM2','NCC','CNRM-CERFACS','IPSL','MPI']\n",
    "rcm_list = ['CLMcom','SMHI-RCA4','KNMI','DMI-HIRHAM5','REMO2015','ALADIN63','REMO2009','WRF381P']"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "id": "d878804e-5b29-422a-8396-a04fdc6da9f6",
   "metadata": {
    "tags": []
   },
   "outputs": [],
   "source": [
    "def mean_bias_tas(string,eobs=eobs_tas,coords=['longitude','latitude']):\n",
    "    model=xr.open_mfdataset(string)\n",
    "    model = model.convert_calendar('noleap',align_on='year')\n",
    "    try:\n",
    "        model = model.reset_coords('height')\n",
    "    except ValueError:\n",
    "        print('no height')\n",
    "    model = model.assign_coords({coord: model[coord].round(3) for coord in coords})\n",
    "    mean_bias = model.tas[0] -273.15 - eobs.tg[0]\n",
    "    return mean_bias"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "id": "c48a60fe-9ce0-4f18-add4-31cf3b49354c",
   "metadata": {},
   "outputs": [],
   "source": [
    "def mean_bias_pr(string,eobs=eobs_pr,coords=['longitude','latitude']):\n",
    "    model=xr.open_mfdataset(string)\n",
    "    model = model.convert_calendar('noleap',align_on='year')\n",
    "    model = model.assign_coords({coord: model[coord].round(3) for coord in coords})\n",
    "    mean_bias = (model.pr[0]*86400 - eobs.rr[0]+0.00001)/(eobs.rr[0]+0.00001) *100\n",
    "    return mean_bias"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 8,
   "id": "8713c856-a02f-4bba-a6d5-1801799614ac",
   "metadata": {
    "tags": []
   },
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "100%|████████████████████████████████████████████████████████████████████████████████████| 8/8 [00:00<00:00, 18.24it/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "No files found for EC-EARTH and ALADIN63, moving to next combination.\n",
      "No files found for EC-EARTH and REMO2009, moving to next combination.\n",
      "No files found for EC-EARTH and WRF381P, moving to next combination.\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "100%|████████████████████████████████████████████████████████████████████████████████████| 8/8 [00:00<00:00, 28.18it/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "No files found for HadGEM2 and ALADIN63, moving to next combination.\n",
      "No files found for HadGEM2 and REMO2009, moving to next combination.\n",
      "No files found for HadGEM2 and WRF381P, moving to next combination.\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "100%|████████████████████████████████████████████████████████████████████████████████████| 8/8 [00:00<00:00, 52.74it/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "No files found for NCC and CLMcom, moving to next combination.\n",
      "No files found for NCC and KNMI, moving to next combination.\n",
      "No files found for NCC and ALADIN63, moving to next combination.\n",
      "No files found for NCC and REMO2009, moving to next combination.\n",
      "No files found for NCC and WRF381P, moving to next combination.\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "100%|████████████████████████████████████████████████████████████████████████████████████| 8/8 [00:00<00:00, 34.90it/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "No files found for CNRM-CERFACS and DMI-HIRHAM5, moving to next combination.\n",
      "No files found for CNRM-CERFACS and REMO2015, moving to next combination.\n",
      "No files found for CNRM-CERFACS and REMO2009, moving to next combination.\n",
      "No files found for CNRM-CERFACS and WRF381P, moving to next combination.\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "100%|████████████████████████████████████████████████████████████████████████████████████| 8/8 [00:00<00:00, 76.71it/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "No files found for IPSL and CLMcom, moving to next combination.\n",
      "No files found for IPSL and KNMI, moving to next combination.\n",
      "No files found for IPSL and DMI-HIRHAM5, moving to next combination.\n",
      "No files found for IPSL and REMO2015, moving to next combination.\n",
      "No files found for IPSL and ALADIN63, moving to next combination.\n",
      "No files found for IPSL and REMO2009, moving to next combination.\n",
      "no height\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "100%|████████████████████████████████████████████████████████████████████████████████████| 8/8 [00:00<00:00, 43.93it/s]"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "No files found for MPI and KNMI, moving to next combination.\n",
      "No files found for MPI and DMI-HIRHAM5, moving to next combination.\n",
      "No files found for MPI and REMO2015, moving to next combination.\n",
      "No files found for MPI and ALADIN63, moving to next combination.\n",
      "No files found for MPI and WRF381P, moving to next combination.\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\n"
     ]
    }
   ],
   "source": [
    "#calculate tas bias in every grid cell\n",
    "bias_list_tas=[]\n",
    "#coords = ['rlon','rlat']\n",
    "\n",
    "names_list = []\n",
    "\n",
    "for gcm in gcm_list:\n",
    "    for rcm in tqdm(rcm_list):\n",
    "        try:\n",
    "            string=f'models_faizan/all_projections/tas*{gcm}*{rcm}*1971_2000_timmean.nc4'\n",
    "            files = glob.glob(string)\n",
    "            bias_list_tas.append(mean_bias_tas(files))\n",
    "            names_list.append(str(gcm)+'_'+str(rcm))\n",
    "        except OSError:\n",
    "            print(f\"No files found for {gcm} and {rcm}, moving to next combination.\")\n",
    "            continue\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "id": "4807f37d-1b8f-4983-81a7-2b17ae6fbaf4",
   "metadata": {
    "tags": []
   },
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "100%|████████████████████████████████████████████████████████████████████████████████████| 8/8 [00:00<00:00, 22.24it/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "No files found for EC-EARTH and ALADIN63, moving to next combination.\n",
      "No files found for EC-EARTH and REMO2009, moving to next combination.\n",
      "No files found for EC-EARTH and WRF381P, moving to next combination.\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "100%|████████████████████████████████████████████████████████████████████████████████████| 8/8 [00:00<00:00, 23.57it/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "No files found for HadGEM2 and ALADIN63, moving to next combination.\n",
      "No files found for HadGEM2 and REMO2009, moving to next combination.\n",
      "No files found for HadGEM2 and WRF381P, moving to next combination.\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      " 50%|██████████████████████████████████████████                                          | 4/8 [00:00<00:00, 37.01it/s]"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "No files found for NCC and CLMcom, moving to next combination.\n",
      "No files found for NCC and KNMI, moving to next combination.\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "100%|████████████████████████████████████████████████████████████████████████████████████| 8/8 [00:00<00:00, 38.01it/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "No files found for NCC and ALADIN63, moving to next combination.\n",
      "No files found for NCC and REMO2009, moving to next combination.\n",
      "No files found for NCC and WRF381P, moving to next combination.\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      " 25%|█████████████████████                                                               | 2/8 [00:00<00:00, 16.94it/s]"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "No files found for CNRM-CERFACS and DMI-HIRHAM5, moving to next combination.\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "100%|████████████████████████████████████████████████████████████████████████████████████| 8/8 [00:00<00:00, 32.37it/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "No files found for CNRM-CERFACS and REMO2015, moving to next combination.\n",
      "No files found for CNRM-CERFACS and REMO2009, moving to next combination.\n",
      "No files found for CNRM-CERFACS and WRF381P, moving to next combination.\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "100%|████████████████████████████████████████████████████████████████████████████████████| 8/8 [00:00<00:00, 71.84it/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "No files found for IPSL and CLMcom, moving to next combination.\n",
      "No files found for IPSL and KNMI, moving to next combination.\n",
      "No files found for IPSL and DMI-HIRHAM5, moving to next combination.\n",
      "No files found for IPSL and REMO2015, moving to next combination.\n",
      "No files found for IPSL and ALADIN63, moving to next combination.\n",
      "No files found for IPSL and REMO2009, moving to next combination.\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "100%|████████████████████████████████████████████████████████████████████████████████████| 8/8 [00:00<00:00, 44.88it/s]\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "No files found for MPI and KNMI, moving to next combination.\n",
      "No files found for MPI and DMI-HIRHAM5, moving to next combination.\n",
      "No files found for MPI and REMO2015, moving to next combination.\n",
      "No files found for MPI and ALADIN63, moving to next combination.\n",
      "No files found for MPI and WRF381P, moving to next combination.\n"
     ]
    }
   ],
   "source": [
    "#calculate pr bias in every grid cell\n",
    "bias_list_pr=[]\n",
    "#coords = ['rlon','rlat']\n",
    "\n",
    "for gcm in gcm_list:\n",
    "    for rcm in tqdm(rcm_list):\n",
    "        try:\n",
    "            string=f'models_faizan/all_projections/pr*{gcm}*{rcm}*1971_2000_timmean.nc4'\n",
    "            bias_list_pr.append(mean_bias_pr(string))    \n",
    "        except OSError:\n",
    "            print(f\"No files found for {gcm} and {rcm}, moving to next combination.\")\n",
    "            continue"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 10,
   "id": "c8db0227-3248-4b7d-9512-a46dc0f1950d",
   "metadata": {
    "collapsed": true,
    "jupyter": {
     "outputs_hidden": true
    },
    "tags": []
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n",
      "Skipping region 141 due to error: No data found in bounds.\n",
      "Skipping region 142 due to error: No data found in bounds.\n",
      "Skipping region 143 due to error: No data found in bounds.\n",
      "Skipping region 144 due to error: No data found in bounds.\n",
      "Skipping region 145 due to error: No data found in bounds.\n"
     ]
    }
   ],
   "source": [
    "#creates lists of bias for every NUTS2 region, skips region outside of CORDEX domain\n",
    "\n",
    "regions_gdf = gpd.read_file('NUTS_2_SHAPE.shp', encoding='utf-8')\n",
    "crs = regions_gdf.crs\n",
    "\n",
    "\n",
    "#crs = ccrs.RotatedPole(pole_latitude=39.25, pole_longitude=-162)\n",
    "\n",
    "#regions_gdf=regions_gdf.to_crs(crs)\n",
    "\n",
    "regions_ex=regions_gdf.explode(column='NUTS_ID')\n",
    "\n",
    "tas_bias_masked_list=[]\n",
    "for model in np.arange(len(bias_list_tas)):\n",
    "    for i in np.arange(len(regions_gdf)):\n",
    "        try:\n",
    "            bias_list_tas[model] = bias_list_tas[model].rio.write_crs(crs)\n",
    "            cropped = bias_list_tas[model].rio.clip(geometries=[regions_ex['geometry'][i]], all_touched=True)\n",
    "            tas_bias_masked_list.append(cropped)\n",
    "        except Exception as e:\n",
    "            print(f\"Skipping region {i} due to error: {e}\")\n",
    "            continue\n",
    "        #bias_median[i,model] = np.median(bias_masked)\n",
    "        #bias_iqr[i,model] = np.quantile(bias_masked,0.75,method='midpoint') - np.quantile(bias_masked,0.25,method='midpoint')\n",
    "        \n",
    "pr_bias_masked_list=[]\n",
    "for model in np.arange(len(bias_list_pr)):\n",
    "    for i in np.arange(len(regions_gdf)):\n",
    "        try:\n",
    "            bias_list_pr[model] = bias_list_pr[model].rio.write_crs(crs)\n",
    "            cropped = bias_list_pr[model].rio.clip(geometries=[regions_ex['geometry'][i]], all_touched=True)\n",
    "            pr_bias_masked_list.append(cropped)\n",
    "        except Exception as e:\n",
    "            print(f\"Skipping region {i} due to error: {e}\")\n",
    "            continue"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 11,
   "id": "820417f1-51b8-49b6-bca3-2f732755ee97",
   "metadata": {
    "tags": []
   },
   "outputs": [],
   "source": [
    "#calculate median bias inside each region\n",
    "\n",
    "weights = [ np.cos(np.deg2rad(i.latitude)) for i in tas_bias_masked_list ]\n",
    "\n",
    "combinations=20 ######## UPDATE this to become automatic i.e., calculate automatically how many gcm-rcm combinations exist\n",
    "\n",
    "\n",
    "median_bias_tas = [i.weighted(weights[num]).quantile(0.5, dim=('latitude', 'longitude')).values for i,num in zip(tas_bias_masked_list,range(len(tas_bias_masked_list)))]\n",
    "#median_bias_tas = np.array(median_bias_tas)\n",
    "#median_bias_tas = median_bias_tas.reshape(combinations,-1)\n",
    "\n",
    "median_bias_pr = [i.weighted(weights[num]).quantile(0.5, dim=('latitude', 'longitude')).values for i,num in zip(pr_bias_masked_list,range(len(pr_bias_masked_list)))]\n",
    "#median_bias_pr = np.array(median_bias_pr)\n",
    "#median_bias_pr = median_bias_pr.reshape(combinations,-1)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "5dd4ff42-dfff-4761-8015-0867e1315ff8",
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "markdown",
   "id": "89e17430-7875-48c3-a020-370811eceb7e",
   "metadata": {},
   "source": [
    "## Create map and attach plots to regions"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 12,
   "id": "3c8515ad-640d-42c6-9946-014d0fef1a5b",
   "metadata": {
    "tags": []
   },
   "outputs": [],
   "source": [
    "regions_gdf = gpd.read_file('NUTS_2_SHAPE.shp', encoding='utf-8')\n",
    "\n",
    "regions_gdf = regions_gdf['geometry']\n",
    "\n",
    "df=gpd.GeoSeries.get_coordinates(regions_gdf)\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 13,
   "id": "6a509c02-47c6-4b39-aabb-3af0aab31429",
   "metadata": {
    "tags": []
   },
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "100%|████████████████████████████████████████████████████████████████████████████████| 329/329 [00:28<00:00, 11.64it/s]\n"
     ]
    }
   ],
   "source": [
    "\n",
    "#drop extra-continental regions\n",
    "#regions_ex2= regions_ex.drop(index=[141,142,143,144,145,163,256,257]).reset_index(drop=True)\n",
    "\n",
    "regions_ex2= regions_ex.drop(index=[141,142,143,144,145]).reset_index(drop=True)\n",
    "\n",
    "#empty list to be filled with all plots produced\n",
    "plots_list = []\n",
    "\n",
    "#label for each model combination\n",
    "label_list = names_list\n",
    "\n",
    "cmap=plt.get_cmap('tab20')\n",
    "color_list = cmap(np.linspace(0, 1,combinations+1))\n",
    "\n",
    "\n",
    "#loop through all regions\n",
    "for region in tqdm(range(len(regions_ex2))):\n",
    "\n",
    "\n",
    "    #make one figure for each region\n",
    "    fig,ax=plt.subplots(1,1,figsize=(6,4))\n",
    "    #region_bias_tas = median_bias_tas[:median_bias_tas/combinations]\n",
    "    #region_bias_pr = median_bias_pr[:,region]\n",
    "    \n",
    "    #in the same figure for each region, plot temp bias vs prec bias, with each color representing one model and having a label\n",
    "    for model in range(combinations):\n",
    "        index = region + ( int( len(median_bias_tas)/combinations) * model)\n",
    "        region_plot=ax.plot(median_bias_tas[index],median_bias_pr[index],'o',color=color_list[model],label=label_list[model])    \n",
    "        \n",
    "    \n",
    "    #Add 0,0 marker for observations\n",
    "    ax.plot(0,0,'s',color='k',markersize=7,label='Observations')\n",
    "    \n",
    "    #Axes titles\n",
    "    ax.set_ylabel('Mean precipitation bias (%)')\n",
    "    ax.set_xlabel('Mean temperature bias (°C)')\n",
    "    \n",
    "    #Set origin to 0,0 and middle of the plot\n",
    "    xlim = ax.get_xlim()\n",
    "    ylim = ax.get_ylim()\n",
    "    ax.plot(xlim, [0, 0], color='b', linewidth=1.5,label='0 Prec bias')\n",
    "    ax.plot([0, 0], ylim, color='r', linewidth=1.5,label='0 Temp bias')\n",
    "    \n",
    "    \n",
    "    #Add legend and append the plot for each region to the list\n",
    "    handles, labels = ax.get_legend_handles_labels()\n",
    "    by_label = dict(zip(labels, handles))\n",
    "    ax.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.05, 0.95), loc='upper left', fontsize=6)\n",
    "    \n",
    "    \n",
    "    plt.title(regions_ex2['NUTS_NAME'][region])\n",
    "    \n",
    "    fig.subplots_adjust(right=0.95)\n",
    "    plt.tight_layout()\n",
    "    \n",
    "    plots_list.append(fig)\n",
    "    plt.close(fig)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 15,
   "id": "9a1e53f4-1047-43f9-859a-31120137064e",
   "metadata": {
    "tags": []
   },
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "100%|████████████████████████████████████████████████████████████████████████████████| 329/329 [00:09<00:00, 33.12it/s]\n"
     ]
    }
   ],
   "source": [
    "#### PLOTLY\n",
    "\n",
    "label_list = names_list\n",
    "\n",
    "# Drop extra-continental regions\n",
    "regions_ex2 = regions_ex.drop(index=[141, 142, 143, 144, 145]).reset_index(drop=True)\n",
    "\n",
    "# Empty list to store plotly figures\n",
    "plots_list_plotly = []\n",
    "\n",
    "# Color list for each model combination\n",
    "colorscale = pcolors.qualitative.Light24\n",
    "color_list = [pcolors.sample_colorscale(colorscale, i / combinations) for i in range(combinations + 1)]\n",
    "\n",
    "# Loop through all regions\n",
    "for region in tqdm(range(len(regions_ex2))):\n",
    "    region_name = regions_ex2['NUTS_NAME'][region]\n",
    "    \n",
    "    # Data collection for plotting\n",
    "    x_data = []\n",
    "    y_data = []\n",
    "    model_labels = []\n",
    "    colors = []\n",
    "\n",
    "    for model in range(combinations):\n",
    "        index = region + (int(len(median_bias_tas) / combinations) * model)\n",
    "        x_data.append(median_bias_tas[index])\n",
    "        y_data.append(median_bias_pr[index])\n",
    "        model_labels.append(label_list[model])\n",
    "        colors.append(color_list[model])\n",
    "\n",
    "    # Create Plotly figure\n",
    "    fig = go.Figure()\n",
    "\n",
    "    # Add model data points\n",
    "    for i in range(combinations):\n",
    "        fig.add_trace(go.Scatter(\n",
    "            x=[x_data[i]],\n",
    "            y=[y_data[i]],\n",
    "            mode='markers',\n",
    "            marker=dict(color=colors[i], size=10),\n",
    "            name=model_labels[i],\n",
    "            customdata=[[region, model_labels[i]]],  # For interaction\n",
    "            hovertemplate=(\n",
    "                f\"<b>Model:</b> {model_labels[i]}<br>\"\n",
    "                f\"<b>Temp Bias:</b> {str(x_data[i])[:4]} °C<br>\"\n",
    "                f\"<b>Prec Bias:</b> {str(y_data[i])[:4]} %<extra></extra>\"\n",
    "            )\n",
    "        ))\n",
    "\n",
    "    # Add origin point (observations)\n",
    "    fig.add_trace(go.Scatter(\n",
    "        x=[0], y=[0],\n",
    "        mode='markers',\n",
    "        marker=dict(color='black', size=10, symbol='square'),\n",
    "        name='Observations',\n",
    "        hoverinfo='skip'\n",
    "    ))\n",
    "\n",
    "    # Add reference lines at zero\n",
    "    fig.add_trace(go.Scatter(x=[min(x_data + [0]), max(x_data + [0])], y=[0, 0],\n",
    "    mode='lines',\n",
    "    line=dict(color='blue', width=1.5),\n",
    "    name='0 Prec bias',\n",
    "    showlegend=True))\n",
    "\n",
    "\n",
    "    fig.add_trace(go.Scatter(x=[0, 0], y=[min(y_data + [0]), max(y_data + [0])],\n",
    "    mode='lines',\n",
    "    line=dict(color='red', width=1.5),\n",
    "    name='0 Temp bias',\n",
    "    showlegend=True))\n",
    "\n",
    "    # Layout\n",
    "    fig.update_layout(\n",
    "        title=f\"{region_name}\",\n",
    "        xaxis_title=\"Mean temperature bias (°C)\",\n",
    "        yaxis_title=\"Mean precipitation bias (%)\",\n",
    "        legend=dict(font=dict(size=8), x=1.02, y=0.95),\n",
    "        margin=dict(l=40, r=40, t=40, b=40),\n",
    "        height=400,\n",
    "        width=600)\n",
    "\n",
    "    \n",
    "    # Append to list\n",
    "    plots_list_plotly.append(fig)   \n",
    "    "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 66,
   "id": "297cc08c-0f16-4e25-8886-925e1a625d49",
   "metadata": {
    "tags": []
   },
   "outputs": [
    {
     "data": {
      "application/vnd.plotly.v1+json": {
       "config": {
        "plotlyServerURL": "https://plot.ly"
       },
       "data": [
        {
         "customdata": [
          [
           100,
           "EC-EARTH_CLMcom"
          ]
         ],
         "hovertemplate": "<b>Model:</b> EC-EARTH_CLMcom<br><b>Temp Bias:</b> -0.7 °C<br><b>Prec Bias:</b> 30.6 %<extra></extra>",
         "marker": {
          "color": [
           "rgb(253, 50, 22)"
          ],
          "size": 10
         },
         "mode": "markers",
         "name": "EC-EARTH_CLMcom",
         "type": "scatter",
         "x": [
          -0.7552823990843343
         ],
         "y": [
          30.676476697870804
         ]
        },
        {
         "customdata": [
          [
           100,
           "EC-EARTH_SMHI-RCA4"
          ]
         ],
         "hovertemplate": "<b>Model:</b> EC-EARTH_SMHI-RCA4<br><b>Temp Bias:</b> -1.4 °C<br><b>Prec Bias:</b> -7.7 %<extra></extra>",
         "marker": {
          "color": [
           "rgb(16, 234, 83)"
          ],
          "size": 10
         },
         "mode": "markers",
         "name": "EC-EARTH_SMHI-RCA4",
         "type": "scatter",
         "x": [
          -1.4171872210163374
         ],
         "y": [
          -7.754017444470064
         ]
        },
        {
         "customdata": [
          [
           100,
           "EC-EARTH_KNMI"
          ]
         ],
         "hovertemplate": "<b>Model:</b> EC-EARTH_KNMI<br><b>Temp Bias:</b> -2.7 °C<br><b>Prec Bias:</b> 20.9 %<extra></extra>",
         "marker": {
          "color": [
           "rgb(150, 146, 235)"
          ],
          "size": 10
         },
         "mode": "markers",
         "name": "EC-EARTH_KNMI",
         "type": "scatter",
         "x": [
          -2.756673228673225
         ],
         "y": [
          20.9153754385099
         ]
        },
        {
         "customdata": [
          [
           100,
           "EC-EARTH_DMI-HIRHAM5"
          ]
         ],
         "hovertemplate": "<b>Model:</b> EC-EARTH_DMI-HIRHAM5<br><b>Temp Bias:</b> -1.2 °C<br><b>Prec Bias:</b> 3.65 %<extra></extra>",
         "marker": {
          "color": [
           "rgb(254, 117, 200)"
          ],
          "size": 10
         },
         "mode": "markers",
         "name": "EC-EARTH_DMI-HIRHAM5",
         "type": "scatter",
         "x": [
          -1.2887597419507424
         ],
         "y": [
          3.654807816569974
         ]
        },
        {
         "customdata": [
          [
           100,
           "EC-EARTH_REMO2015"
          ]
         ],
         "hovertemplate": "<b>Model:</b> EC-EARTH_REMO2015<br><b>Temp Bias:</b> -0.1 °C<br><b>Prec Bias:</b> 2.19 %<extra></extra>",
         "marker": {
          "color": [
           "rgb(109, 149, 235)"
          ],
          "size": 10
         },
         "mode": "markers",
         "name": "EC-EARTH_REMO2015",
         "type": "scatter",
         "x": [
          -0.1417296636306438
         ],
         "y": [
          2.1957434757533427
         ]
        },
        {
         "customdata": [
          [
           100,
           "HadGEM2_CLMcom"
          ]
         ],
         "hovertemplate": "<b>Model:</b> HadGEM2_CLMcom<br><b>Temp Bias:</b> 0.76 °C<br><b>Prec Bias:</b> -3.0 %<extra></extra>",
         "marker": {
          "color": [
           "rgb(188, 249, 92)"
          ],
          "size": 10
         },
         "mode": "markers",
         "name": "HadGEM2_CLMcom",
         "type": "scatter",
         "x": [
          0.7647336183236615
         ],
         "y": [
          -3.09269158075561
         ]
        },
        {
         "customdata": [
          [
           100,
           "HadGEM2_SMHI-RCA4"
          ]
         ],
         "hovertemplate": "<b>Model:</b> HadGEM2_SMHI-RCA4<br><b>Temp Bias:</b> 0.25 °C<br><b>Prec Bias:</b> -14. %<extra></extra>",
         "marker": {
          "color": [
           "rgb(254, 160, 24)"
          ],
          "size": 10
         },
         "mode": "markers",
         "name": "HadGEM2_SMHI-RCA4",
         "type": "scatter",
         "x": [
          0.25594916290833236
         ],
         "y": [
          -14.556739887662903
         ]
        },
        {
         "customdata": [
          [
           100,
           "HadGEM2_KNMI"
          ]
         ],
         "hovertemplate": "<b>Model:</b> HadGEM2_KNMI<br><b>Temp Bias:</b> -1.1 °C<br><b>Prec Bias:</b> 18.3 %<extra></extra>",
         "marker": {
          "color": [
           "rgb(79, 156, 93)"
          ],
          "size": 10
         },
         "mode": "markers",
         "name": "HadGEM2_KNMI",
         "type": "scatter",
         "x": [
          -1.1713855548760241
         ],
         "y": [
          18.301829327955005
         ]
        },
        {
         "customdata": [
          [
           100,
           "HadGEM2_DMI-HIRHAM5"
          ]
         ],
         "hovertemplate": "<b>Model:</b> HadGEM2_DMI-HIRHAM5<br><b>Temp Bias:</b> 0.60 °C<br><b>Prec Bias:</b> -10. %<extra></extra>",
         "marker": {
          "color": [
           "rgb(234, 150, 226)"
          ],
          "size": 10
         },
         "mode": "markers",
         "name": "HadGEM2_DMI-HIRHAM5",
         "type": "scatter",
         "x": [
          0.6049662037405695
         ],
         "y": [
          -10.057160748952777
         ]
        },
        {
         "customdata": [
          [
           100,
           "HadGEM2_REMO2015"
          ]
         ],
         "hovertemplate": "<b>Model:</b> HadGEM2_REMO2015<br><b>Temp Bias:</b> 1.55 °C<br><b>Prec Bias:</b> -1.4 %<extra></extra>",
         "marker": {
          "color": [
           "rgb(218, 70, 171)"
          ],
          "size": 10
         },
         "mode": "markers",
         "name": "HadGEM2_REMO2015",
         "type": "scatter",
         "x": [
          1.5533013323294933
         ],
         "y": [
          -1.4203766768350035
         ]
        },
        {
         "customdata": [
          [
           100,
           "NCC_SMHI-RCA4"
          ]
         ],
         "hovertemplate": "<b>Model:</b> NCC_SMHI-RCA4<br><b>Temp Bias:</b> 0.09 °C<br><b>Prec Bias:</b> -3.3 %<extra></extra>",
         "marker": {
          "color": [
           "rgb(162, 88, 206)"
          ],
          "size": 10
         },
         "mode": "markers",
         "name": "NCC_SMHI-RCA4",
         "type": "scatter",
         "x": [
          0.09219400808735957
         ],
         "y": [
          -3.331383563240392
         ]
        },
        {
         "customdata": [
          [
           100,
           "NCC_DMI-HIRHAM5"
          ]
         ],
         "hovertemplate": "<b>Model:</b> NCC_DMI-HIRHAM5<br><b>Temp Bias:</b> 0.47 °C<br><b>Prec Bias:</b> -3.7 %<extra></extra>",
         "marker": {
          "color": [
           "rgb(38, 166, 215)"
          ],
          "size": 10
         },
         "mode": "markers",
         "name": "NCC_DMI-HIRHAM5",
         "type": "scatter",
         "x": [
          0.4794859435892186
         ],
         "y": [
          -3.769711923097435
         ]
        },
        {
         "customdata": [
          [
           100,
           "NCC_REMO2015"
          ]
         ],
         "hovertemplate": "<b>Model:</b> NCC_REMO2015<br><b>Temp Bias:</b> 0.90 °C<br><b>Prec Bias:</b> -3.9 %<extra></extra>",
         "marker": {
          "color": [
           "rgb(146, 150, 49)"
          ],
          "size": 10
         },
         "mode": "markers",
         "name": "NCC_REMO2015",
         "type": "scatter",
         "x": [
          0.9010910335979359
         ],
         "y": [
          -3.9392930754865247
         ]
        },
        {
         "customdata": [
          [
           100,
           "CNRM-CERFACS_KNMI"
          ]
         ],
         "hovertemplate": "<b>Model:</b> CNRM-CERFACS_KNMI<br><b>Temp Bias:</b> -2.6 °C<br><b>Prec Bias:</b> 50.8 %<extra></extra>",
         "marker": {
          "color": [
           "rgb(200, 246, 218)"
          ],
          "size": 10
         },
         "mode": "markers",
         "name": "CNRM-CERFACS_KNMI",
         "type": "scatter",
         "x": [
          -2.621824679647758
         ],
         "y": [
          50.86273542936121
         ]
        },
        {
         "customdata": [
          [
           100,
           "CNRM-CERFACS_ALADIN63"
          ]
         ],
         "hovertemplate": "<b>Model:</b> CNRM-CERFACS_ALADIN63<br><b>Temp Bias:</b> -1.7 °C<br><b>Prec Bias:</b> 70.3 %<extra></extra>",
         "marker": {
          "color": [
           "rgb(233, 25, 148)"
          ],
          "size": 10
         },
         "mode": "markers",
         "name": "CNRM-CERFACS_ALADIN63",
         "type": "scatter",
         "x": [
          -1.7337617843174475
         ],
         "y": [
          70.39465715222727
         ]
        },
        {
         "customdata": [
          [
           100,
           "IPSL_SMHI-RCA4"
          ]
         ],
         "hovertemplate": "<b>Model:</b> IPSL_SMHI-RCA4<br><b>Temp Bias:</b> -0.4 °C<br><b>Prec Bias:</b> 6.12 %<extra></extra>",
         "marker": {
          "color": [
           "rgb(82, 251, 165)"
          ],
          "size": 10
         },
         "mode": "markers",
         "name": "IPSL_SMHI-RCA4",
         "type": "scatter",
         "x": [
          -0.46861886151068016
         ],
         "y": [
          6.128590937718343
         ]
        },
        {
         "customdata": [
          [
           100,
           "IPSL_WRF381P"
          ]
         ],
         "hovertemplate": "<b>Model:</b> IPSL_WRF381P<br><b>Temp Bias:</b> -0.8 °C<br><b>Prec Bias:</b> 32.7 %<extra></extra>",
         "marker": {
          "color": [
           "rgb(190, 225, 95)"
          ],
          "size": 10
         },
         "mode": "markers",
         "name": "IPSL_WRF381P",
         "type": "scatter",
         "x": [
          -0.8971256678551027
         ],
         "y": [
          32.75942321691502
         ]
        },
        {
         "customdata": [
          [
           100,
           "MPI_CLMcom"
          ]
         ],
         "hovertemplate": "<b>Model:</b> MPI_CLMcom<br><b>Temp Bias:</b> -0.0 °C<br><b>Prec Bias:</b> 54.3 %<extra></extra>",
         "marker": {
          "color": [
           "rgb(164, 155, 82)"
          ],
          "size": 10
         },
         "mode": "markers",
         "name": "MPI_CLMcom",
         "type": "scatter",
         "x": [
          -0.012765430399688908
         ],
         "y": [
          54.34185803902375
         ]
        },
        {
         "customdata": [
          [
           100,
           "MPI_SMHI-RCA4"
          ]
         ],
         "hovertemplate": "<b>Model:</b> MPI_SMHI-RCA4<br><b>Temp Bias:</b> 0.04 °C<br><b>Prec Bias:</b> 16.4 %<extra></extra>",
         "marker": {
          "color": [
           "rgb(145, 121, 189)"
          ],
          "size": 10
         },
         "mode": "markers",
         "name": "MPI_SMHI-RCA4",
         "type": "scatter",
         "x": [
          0.04996812031626846
         ],
         "y": [
          16.45607454101971
         ]
        },
        {
         "customdata": [
          [
           100,
           "MPI_REMO2009"
          ]
         ],
         "hovertemplate": "<b>Model:</b> MPI_REMO2009<br><b>Temp Bias:</b> 0.85 °C<br><b>Prec Bias:</b> 4.23 %<extra></extra>",
         "marker": {
          "color": [
           "rgb(233, 108, 103)"
          ],
          "size": 10
         },
         "mode": "markers",
         "name": "MPI_REMO2009",
         "type": "scatter",
         "x": [
          0.8586443612603661
         ],
         "y": [
          4.236100432615828
         ]
        },
        {
         "hoverinfo": "skip",
         "marker": {
          "color": "black",
          "size": 10,
          "symbol": "square"
         },
         "mode": "markers",
         "name": "Observations",
         "type": "scatter",
         "x": [
          0
         ],
         "y": [
          0
         ]
        },
        {
         "line": {
          "color": "blue",
          "width": 1.5
         },
         "mode": "lines",
         "name": "0 Prec bias",
         "showlegend": true,
         "type": "scatter",
         "x": [
          -2.756673228673225,
          1.5533013323294933
         ],
         "y": [
          0,
          0
         ]
        },
        {
         "line": {
          "color": "red",
          "width": 1.5
         },
         "mode": "lines",
         "name": "0 Temp bias",
         "showlegend": true,
         "type": "scatter",
         "x": [
          0,
          0
         ],
         "y": [
          -14.556739887662903,
          70.39465715222727
         ]
        }
       ],
       "layout": {
        "height": 400,
        "legend": {
         "font": {
          "size": 8
         },
         "x": 1.02,
         "y": 0.95
        },
        "margin": {
         "b": 40,
         "l": 40,
         "r": 40,
         "t": 40
        },
        "template": {
         "data": {
          "bar": [
           {
            "error_x": {
             "color": "#2a3f5f"
            },
            "error_y": {
             "color": "#2a3f5f"
            },
            "marker": {
             "line": {
              "color": "#E5ECF6",
              "width": 0.5
             },
             "pattern": {
              "fillmode": "overlay",
              "size": 10,
              "solidity": 0.2
             }
            },
            "type": "bar"
           }
          ],
          "barpolar": [
           {
            "marker": {
             "line": {
              "color": "#E5ECF6",
              "width": 0.5
             },
             "pattern": {
              "fillmode": "overlay",
              "size": 10,
              "solidity": 0.2
             }
            },
            "type": "barpolar"
           }
          ],
          "carpet": [
           {
            "aaxis": {
             "endlinecolor": "#2a3f5f",
             "gridcolor": "white",
             "linecolor": "white",
             "minorgridcolor": "white",
             "startlinecolor": "#2a3f5f"
            },
            "baxis": {
             "endlinecolor": "#2a3f5f",
             "gridcolor": "white",
             "linecolor": "white",
             "minorgridcolor": "white",
             "startlinecolor": "#2a3f5f"
            },
            "type": "carpet"
           }
          ],
          "choropleth": [
           {
            "colorbar": {
             "outlinewidth": 0,
             "ticks": ""
            },
            "type": "choropleth"
           }
          ],
          "contour": [
           {
            "colorbar": {
             "outlinewidth": 0,
             "ticks": ""
            },
            "colorscale": [
             [
              0,
              "#0d0887"
             ],
             [
              0.1111111111111111,
              "#46039f"
             ],
             [
              0.2222222222222222,
              "#7201a8"
             ],
             [
              0.3333333333333333,
              "#9c179e"
             ],
             [
              0.4444444444444444,
              "#bd3786"
             ],
             [
              0.5555555555555556,
              "#d8576b"
             ],
             [
              0.6666666666666666,
              "#ed7953"
             ],
             [
              0.7777777777777778,
              "#fb9f3a"
             ],
             [
              0.8888888888888888,
              "#fdca26"
             ],
             [
              1,
              "#f0f921"
             ]
            ],
            "type": "contour"
           }
          ],
          "contourcarpet": [
           {
            "colorbar": {
             "outlinewidth": 0,
             "ticks": ""
            },
            "type": "contourcarpet"
           }
          ],
          "heatmap": [
           {
            "colorbar": {
             "outlinewidth": 0,
             "ticks": ""
            },
            "colorscale": [
             [
              0,
              "#0d0887"
             ],
             [
              0.1111111111111111,
              "#46039f"
             ],
             [
              0.2222222222222222,
              "#7201a8"
             ],
             [
              0.3333333333333333,
              "#9c179e"
             ],
             [
              0.4444444444444444,
              "#bd3786"
             ],
             [
              0.5555555555555556,
              "#d8576b"
             ],
             [
              0.6666666666666666,
              "#ed7953"
             ],
             [
              0.7777777777777778,
              "#fb9f3a"
             ],
             [
              0.8888888888888888,
              "#fdca26"
             ],
             [
              1,
              "#f0f921"
             ]
            ],
            "type": "heatmap"
           }
          ],
          "heatmapgl": [
           {
            "colorbar": {
             "outlinewidth": 0,
             "ticks": ""
            },
            "colorscale": [
             [
              0,
              "#0d0887"
             ],
             [
              0.1111111111111111,
              "#46039f"
             ],
             [
              0.2222222222222222,
              "#7201a8"
             ],
             [
              0.3333333333333333,
              "#9c179e"
             ],
             [
              0.4444444444444444,
              "#bd3786"
             ],
             [
              0.5555555555555556,
              "#d8576b"
             ],
             [
              0.6666666666666666,
              "#ed7953"
             ],
             [
              0.7777777777777778,
              "#fb9f3a"
             ],
             [
              0.8888888888888888,
              "#fdca26"
             ],
             [
              1,
              "#f0f921"
             ]
            ],
            "type": "heatmapgl"
           }
          ],
          "histogram": [
           {
            "marker": {
             "pattern": {
              "fillmode": "overlay",
              "size": 10,
              "solidity": 0.2
             }
            },
            "type": "histogram"
           }
          ],
          "histogram2d": [
           {
            "colorbar": {
             "outlinewidth": 0,
             "ticks": ""
            },
            "colorscale": [
             [
              0,
              "#0d0887"
             ],
             [
              0.1111111111111111,
              "#46039f"
             ],
             [
              0.2222222222222222,
              "#7201a8"
             ],
             [
              0.3333333333333333,
              "#9c179e"
             ],
             [
              0.4444444444444444,
              "#bd3786"
             ],
             [
              0.5555555555555556,
              "#d8576b"
             ],
             [
              0.6666666666666666,
              "#ed7953"
             ],
             [
              0.7777777777777778,
              "#fb9f3a"
             ],
             [
              0.8888888888888888,
              "#fdca26"
             ],
             [
              1,
              "#f0f921"
             ]
            ],
            "type": "histogram2d"
           }
          ],
          "histogram2dcontour": [
           {
            "colorbar": {
             "outlinewidth": 0,
             "ticks": ""
            },
            "colorscale": [
             [
              0,
              "#0d0887"
             ],
             [
              0.1111111111111111,
              "#46039f"
             ],
             [
              0.2222222222222222,
              "#7201a8"
             ],
             [
              0.3333333333333333,
              "#9c179e"
             ],
             [
              0.4444444444444444,
              "#bd3786"
             ],
             [
              0.5555555555555556,
              "#d8576b"
             ],
             [
              0.6666666666666666,
              "#ed7953"
             ],
             [
              0.7777777777777778,
              "#fb9f3a"
             ],
             [
              0.8888888888888888,
              "#fdca26"
             ],
             [
              1,
              "#f0f921"
             ]
            ],
            "type": "histogram2dcontour"
           }
          ],
          "mesh3d": [
           {
            "colorbar": {
             "outlinewidth": 0,
             "ticks": ""
            },
            "type": "mesh3d"
           }
          ],
          "parcoords": [
           {
            "line": {
             "colorbar": {
              "outlinewidth": 0,
              "ticks": ""
             }
            },
            "type": "parcoords"
           }
          ],
          "pie": [
           {
            "automargin": true,
            "type": "pie"
           }
          ],
          "scatter": [
           {
            "fillpattern": {
             "fillmode": "overlay",
             "size": 10,
             "solidity": 0.2
            },
            "type": "scatter"
           }
          ],
          "scatter3d": [
           {
            "line": {
             "colorbar": {
              "outlinewidth": 0,
              "ticks": ""
             }
            },
            "marker": {
             "colorbar": {
              "outlinewidth": 0,
              "ticks": ""
             }
            },
            "type": "scatter3d"
           }
          ],
          "scattercarpet": [
           {
            "marker": {
             "colorbar": {
              "outlinewidth": 0,
              "ticks": ""
             }
            },
            "type": "scattercarpet"
           }
          ],
          "scattergeo": [
           {
            "marker": {
             "colorbar": {
              "outlinewidth": 0,
              "ticks": ""
             }
            },
            "type": "scattergeo"
           }
          ],
          "scattergl": [
           {
            "marker": {
             "colorbar": {
              "outlinewidth": 0,
              "ticks": ""
             }
            },
            "type": "scattergl"
           }
          ],
          "scattermapbox": [
           {
            "marker": {
             "colorbar": {
              "outlinewidth": 0,
              "ticks": ""
             }
            },
            "type": "scattermapbox"
           }
          ],
          "scatterpolar": [
           {
            "marker": {
             "colorbar": {
              "outlinewidth": 0,
              "ticks": ""
             }
            },
            "type": "scatterpolar"
           }
          ],
          "scatterpolargl": [
           {
            "marker": {
             "colorbar": {
              "outlinewidth": 0,
              "ticks": ""
             }
            },
            "type": "scatterpolargl"
           }
          ],
          "scatterternary": [
           {
            "marker": {
             "colorbar": {
              "outlinewidth": 0,
              "ticks": ""
             }
            },
            "type": "scatterternary"
           }
          ],
          "surface": [
           {
            "colorbar": {
             "outlinewidth": 0,
             "ticks": ""
            },
            "colorscale": [
             [
              0,
              "#0d0887"
             ],
             [
              0.1111111111111111,
              "#46039f"
             ],
             [
              0.2222222222222222,
              "#7201a8"
             ],
             [
              0.3333333333333333,
              "#9c179e"
             ],
             [
              0.4444444444444444,
              "#bd3786"
             ],
             [
              0.5555555555555556,
              "#d8576b"
             ],
             [
              0.6666666666666666,
              "#ed7953"
             ],
             [
              0.7777777777777778,
              "#fb9f3a"
             ],
             [
              0.8888888888888888,
              "#fdca26"
             ],
             [
              1,
              "#f0f921"
             ]
            ],
            "type": "surface"
           }
          ],
          "table": [
           {
            "cells": {
             "fill": {
              "color": "#EBF0F8"
             },
             "line": {
              "color": "white"
             }
            },
            "header": {
             "fill": {
              "color": "#C8D4E3"
             },
             "line": {
              "color": "white"
             }
            },
            "type": "table"
           }
          ]
         },
         "layout": {
          "annotationdefaults": {
           "arrowcolor": "#2a3f5f",
           "arrowhead": 0,
           "arrowwidth": 1
          },
          "autotypenumbers": "strict",
          "coloraxis": {
           "colorbar": {
            "outlinewidth": 0,
            "ticks": ""
           }
          },
          "colorscale": {
           "diverging": [
            [
             0,
             "#8e0152"
            ],
            [
             0.1,
             "#c51b7d"
            ],
            [
             0.2,
             "#de77ae"
            ],
            [
             0.3,
             "#f1b6da"
            ],
            [
             0.4,
             "#fde0ef"
            ],
            [
             0.5,
             "#f7f7f7"
            ],
            [
             0.6,
             "#e6f5d0"
            ],
            [
             0.7,
             "#b8e186"
            ],
            [
             0.8,
             "#7fbc41"
            ],
            [
             0.9,
             "#4d9221"
            ],
            [
             1,
             "#276419"
            ]
           ],
           "sequential": [
            [
             0,
             "#0d0887"
            ],
            [
             0.1111111111111111,
             "#46039f"
            ],
            [
             0.2222222222222222,
             "#7201a8"
            ],
            [
             0.3333333333333333,
             "#9c179e"
            ],
            [
             0.4444444444444444,
             "#bd3786"
            ],
            [
             0.5555555555555556,
             "#d8576b"
            ],
            [
             0.6666666666666666,
             "#ed7953"
            ],
            [
             0.7777777777777778,
             "#fb9f3a"
            ],
            [
             0.8888888888888888,
             "#fdca26"
            ],
            [
             1,
             "#f0f921"
            ]
           ],
           "sequentialminus": [
            [
             0,
             "#0d0887"
            ],
            [
             0.1111111111111111,
             "#46039f"
            ],
            [
             0.2222222222222222,
             "#7201a8"
            ],
            [
             0.3333333333333333,
             "#9c179e"
            ],
            [
             0.4444444444444444,
             "#bd3786"
            ],
            [
             0.5555555555555556,
             "#d8576b"
            ],
            [
             0.6666666666666666,
             "#ed7953"
            ],
            [
             0.7777777777777778,
             "#fb9f3a"
            ],
            [
             0.8888888888888888,
             "#fdca26"
            ],
            [
             1,
             "#f0f921"
            ]
           ]
          },
          "colorway": [
           "#636efa",
           "#EF553B",
           "#00cc96",
           "#ab63fa",
           "#FFA15A",
           "#19d3f3",
           "#FF6692",
           "#B6E880",
           "#FF97FF",
           "#FECB52"
          ],
          "font": {
           "color": "#2a3f5f"
          },
          "geo": {
           "bgcolor": "white",
           "lakecolor": "white",
           "landcolor": "#E5ECF6",
           "showlakes": true,
           "showland": true,
           "subunitcolor": "white"
          },
          "hoverlabel": {
           "align": "left"
          },
          "hovermode": "closest",
          "mapbox": {
           "style": "light"
          },
          "paper_bgcolor": "white",
          "plot_bgcolor": "#E5ECF6",
          "polar": {
           "angularaxis": {
            "gridcolor": "white",
            "linecolor": "white",
            "ticks": ""
           },
           "bgcolor": "#E5ECF6",
           "radialaxis": {
            "gridcolor": "white",
            "linecolor": "white",
            "ticks": ""
           }
          },
          "scene": {
           "xaxis": {
            "backgroundcolor": "#E5ECF6",
            "gridcolor": "white",
            "gridwidth": 2,
            "linecolor": "white",
            "showbackground": true,
            "ticks": "",
            "zerolinecolor": "white"
           },
           "yaxis": {
            "backgroundcolor": "#E5ECF6",
            "gridcolor": "white",
            "gridwidth": 2,
            "linecolor": "white",
            "showbackground": true,
            "ticks": "",
            "zerolinecolor": "white"
           },
           "zaxis": {
            "backgroundcolor": "#E5ECF6",
            "gridcolor": "white",
            "gridwidth": 2,
            "linecolor": "white",
            "showbackground": true,
            "ticks": "",
            "zerolinecolor": "white"
           }
          },
          "shapedefaults": {
           "line": {
            "color": "#2a3f5f"
           }
          },
          "ternary": {
           "aaxis": {
            "gridcolor": "white",
            "linecolor": "white",
            "ticks": ""
           },
           "baxis": {
            "gridcolor": "white",
            "linecolor": "white",
            "ticks": ""
           },
           "bgcolor": "#E5ECF6",
           "caxis": {
            "gridcolor": "white",
            "linecolor": "white",
            "ticks": ""
           }
          },
          "title": {
           "x": 0.05
          },
          "xaxis": {
           "automargin": true,
           "gridcolor": "white",
           "linecolor": "white",
           "ticks": "",
           "title": {
            "standoff": 15
           },
           "zerolinecolor": "white",
           "zerolinewidth": 2
          },
          "yaxis": {
           "automargin": true,
           "gridcolor": "white",
           "linecolor": "white",
           "ticks": "",
           "title": {
            "standoff": 15
           },
           "zerolinecolor": "white",
           "zerolinewidth": 2
          }
         }
        },
        "title": {
         "text": "Friuli-Venezia Giulia"
        },
        "width": 600,
        "xaxis": {
         "autorange": true,
         "range": [
          -3.091648908725623,
          1.888277012381891
         ],
         "title": {
          "text": "Mean temperature bias (°C)"
         },
         "type": "linear"
        },
        "yaxis": {
         "autorange": true,
         "range": [
          -21.50497466633167,
          77.34289193089603
         ],
         "title": {
          "text": "Mean precipitation bias (%)"
         },
         "type": "linear"
        }
       }
      },
      "image/png": "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",
      "text/html": [
       "<div>                            <div id=\"58ab1288-036b-4531-8e77-fcbfba46fbde\" class=\"plotly-graph-div\" style=\"height:400px; width:600px;\"></div>            <script type=\"text/javascript\">                require([\"plotly\"], function(Plotly) {                    window.PLOTLYENV=window.PLOTLYENV || {};                                    if (document.getElementById(\"58ab1288-036b-4531-8e77-fcbfba46fbde\")) {                    Plotly.newPlot(                        \"58ab1288-036b-4531-8e77-fcbfba46fbde\",                        [{\"customdata\":[[100,\"EC-EARTH_CLMcom\"]],\"hovertemplate\":\"\\u003cb\\u003eModel:\\u003c\\u002fb\\u003e EC-EARTH_CLMcom\\u003cbr\\u003e\\u003cb\\u003eTemp Bias:\\u003c\\u002fb\\u003e -0.7 \\u00b0C\\u003cbr\\u003e\\u003cb\\u003ePrec Bias:\\u003c\\u002fb\\u003e 30.6 %\\u003cextra\\u003e\\u003c\\u002fextra\\u003e\",\"marker\":{\"color\":[\"rgb(253, 50, 22)\"],\"size\":10},\"mode\":\"markers\",\"name\":\"EC-EARTH_CLMcom\",\"x\":[-0.7552823990843343],\"y\":[30.676476697870804],\"type\":\"scatter\"},{\"customdata\":[[100,\"EC-EARTH_SMHI-RCA4\"]],\"hovertemplate\":\"\\u003cb\\u003eModel:\\u003c\\u002fb\\u003e EC-EARTH_SMHI-RCA4\\u003cbr\\u003e\\u003cb\\u003eTemp Bias:\\u003c\\u002fb\\u003e -1.4 \\u00b0C\\u003cbr\\u003e\\u003cb\\u003ePrec Bias:\\u003c\\u002fb\\u003e -7.7 %\\u003cextra\\u003e\\u003c\\u002fextra\\u003e\",\"marker\":{\"color\":[\"rgb(16, 234, 83)\"],\"size\":10},\"mode\":\"markers\",\"name\":\"EC-EARTH_SMHI-RCA4\",\"x\":[-1.4171872210163374],\"y\":[-7.754017444470064],\"type\":\"scatter\"},{\"customdata\":[[100,\"EC-EARTH_KNMI\"]],\"hovertemplate\":\"\\u003cb\\u003eModel:\\u003c\\u002fb\\u003e EC-EARTH_KNMI\\u003cbr\\u003e\\u003cb\\u003eTemp Bias:\\u003c\\u002fb\\u003e -2.7 \\u00b0C\\u003cbr\\u003e\\u003cb\\u003ePrec Bias:\\u003c\\u002fb\\u003e 20.9 %\\u003cextra\\u003e\\u003c\\u002fextra\\u003e\",\"marker\":{\"color\":[\"rgb(150, 146, 235)\"],\"size\":10},\"mode\":\"markers\",\"name\":\"EC-EARTH_KNMI\",\"x\":[-2.756673228673225],\"y\":[20.9153754385099],\"type\":\"scatter\"},{\"customdata\":[[100,\"EC-EARTH_DMI-HIRHAM5\"]],\"hovertemplate\":\"\\u003cb\\u003eModel:\\u003c\\u002fb\\u003e EC-EARTH_DMI-HIRHAM5\\u003cbr\\u003e\\u003cb\\u003eTemp Bias:\\u003c\\u002fb\\u003e -1.2 \\u00b0C\\u003cbr\\u003e\\u003cb\\u003ePrec Bias:\\u003c\\u002fb\\u003e 3.65 %\\u003cextra\\u003e\\u003c\\u002fextra\\u003e\",\"marker\":{\"color\":[\"rgb(254, 117, 200)\"],\"size\":10},\"mode\":\"markers\",\"name\":\"EC-EARTH_DMI-HIRHAM5\",\"x\":[-1.2887597419507424],\"y\":[3.654807816569974],\"type\":\"scatter\"},{\"customdata\":[[100,\"EC-EARTH_REMO2015\"]],\"hovertemplate\":\"\\u003cb\\u003eModel:\\u003c\\u002fb\\u003e EC-EARTH_REMO2015\\u003cbr\\u003e\\u003cb\\u003eTemp Bias:\\u003c\\u002fb\\u003e -0.1 \\u00b0C\\u003cbr\\u003e\\u003cb\\u003ePrec Bias:\\u003c\\u002fb\\u003e 2.19 %\\u003cextra\\u003e\\u003c\\u002fextra\\u003e\",\"marker\":{\"color\":[\"rgb(109, 149, 235)\"],\"size\":10},\"mode\":\"markers\",\"name\":\"EC-EARTH_REMO2015\",\"x\":[-0.1417296636306438],\"y\":[2.1957434757533427],\"type\":\"scatter\"},{\"customdata\":[[100,\"HadGEM2_CLMcom\"]],\"hovertemplate\":\"\\u003cb\\u003eModel:\\u003c\\u002fb\\u003e HadGEM2_CLMcom\\u003cbr\\u003e\\u003cb\\u003eTemp Bias:\\u003c\\u002fb\\u003e 0.76 \\u00b0C\\u003cbr\\u003e\\u003cb\\u003ePrec Bias:\\u003c\\u002fb\\u003e -3.0 %\\u003cextra\\u003e\\u003c\\u002fextra\\u003e\",\"marker\":{\"color\":[\"rgb(188, 249, 92)\"],\"size\":10},\"mode\":\"markers\",\"name\":\"HadGEM2_CLMcom\",\"x\":[0.7647336183236615],\"y\":[-3.09269158075561],\"type\":\"scatter\"},{\"customdata\":[[100,\"HadGEM2_SMHI-RCA4\"]],\"hovertemplate\":\"\\u003cb\\u003eModel:\\u003c\\u002fb\\u003e HadGEM2_SMHI-RCA4\\u003cbr\\u003e\\u003cb\\u003eTemp Bias:\\u003c\\u002fb\\u003e 0.25 \\u00b0C\\u003cbr\\u003e\\u003cb\\u003ePrec Bias:\\u003c\\u002fb\\u003e -14. %\\u003cextra\\u003e\\u003c\\u002fextra\\u003e\",\"marker\":{\"color\":[\"rgb(254, 160, 24)\"],\"size\":10},\"mode\":\"markers\",\"name\":\"HadGEM2_SMHI-RCA4\",\"x\":[0.25594916290833236],\"y\":[-14.556739887662903],\"type\":\"scatter\"},{\"customdata\":[[100,\"HadGEM2_KNMI\"]],\"hovertemplate\":\"\\u003cb\\u003eModel:\\u003c\\u002fb\\u003e HadGEM2_KNMI\\u003cbr\\u003e\\u003cb\\u003eTemp Bias:\\u003c\\u002fb\\u003e -1.1 \\u00b0C\\u003cbr\\u003e\\u003cb\\u003ePrec Bias:\\u003c\\u002fb\\u003e 18.3 %\\u003cextra\\u003e\\u003c\\u002fextra\\u003e\",\"marker\":{\"color\":[\"rgb(79, 156, 93)\"],\"size\":10},\"mode\":\"markers\",\"name\":\"HadGEM2_KNMI\",\"x\":[-1.1713855548760241],\"y\":[18.301829327955005],\"type\":\"scatter\"},{\"customdata\":[[100,\"HadGEM2_DMI-HIRHAM5\"]],\"hovertemplate\":\"\\u003cb\\u003eModel:\\u003c\\u002fb\\u003e HadGEM2_DMI-HIRHAM5\\u003cbr\\u003e\\u003cb\\u003eTemp Bias:\\u003c\\u002fb\\u003e 0.60 \\u00b0C\\u003cbr\\u003e\\u003cb\\u003ePrec Bias:\\u003c\\u002fb\\u003e -10. %\\u003cextra\\u003e\\u003c\\u002fextra\\u003e\",\"marker\":{\"color\":[\"rgb(234, 150, 226)\"],\"size\":10},\"mode\":\"markers\",\"name\":\"HadGEM2_DMI-HIRHAM5\",\"x\":[0.6049662037405695],\"y\":[-10.057160748952777],\"type\":\"scatter\"},{\"customdata\":[[100,\"HadGEM2_REMO2015\"]],\"hovertemplate\":\"\\u003cb\\u003eModel:\\u003c\\u002fb\\u003e HadGEM2_REMO2015\\u003cbr\\u003e\\u003cb\\u003eTemp Bias:\\u003c\\u002fb\\u003e 1.55 \\u00b0C\\u003cbr\\u003e\\u003cb\\u003ePrec Bias:\\u003c\\u002fb\\u003e -1.4 %\\u003cextra\\u003e\\u003c\\u002fextra\\u003e\",\"marker\":{\"color\":[\"rgb(218, 70, 171)\"],\"size\":10},\"mode\":\"markers\",\"name\":\"HadGEM2_REMO2015\",\"x\":[1.5533013323294933],\"y\":[-1.4203766768350035],\"type\":\"scatter\"},{\"customdata\":[[100,\"NCC_SMHI-RCA4\"]],\"hovertemplate\":\"\\u003cb\\u003eModel:\\u003c\\u002fb\\u003e NCC_SMHI-RCA4\\u003cbr\\u003e\\u003cb\\u003eTemp Bias:\\u003c\\u002fb\\u003e 0.09 \\u00b0C\\u003cbr\\u003e\\u003cb\\u003ePrec Bias:\\u003c\\u002fb\\u003e -3.3 %\\u003cextra\\u003e\\u003c\\u002fextra\\u003e\",\"marker\":{\"color\":[\"rgb(162, 88, 206)\"],\"size\":10},\"mode\":\"markers\",\"name\":\"NCC_SMHI-RCA4\",\"x\":[0.09219400808735957],\"y\":[-3.331383563240392],\"type\":\"scatter\"},{\"customdata\":[[100,\"NCC_DMI-HIRHAM5\"]],\"hovertemplate\":\"\\u003cb\\u003eModel:\\u003c\\u002fb\\u003e NCC_DMI-HIRHAM5\\u003cbr\\u003e\\u003cb\\u003eTemp Bias:\\u003c\\u002fb\\u003e 0.47 \\u00b0C\\u003cbr\\u003e\\u003cb\\u003ePrec Bias:\\u003c\\u002fb\\u003e -3.7 %\\u003cextra\\u003e\\u003c\\u002fextra\\u003e\",\"marker\":{\"color\":[\"rgb(38, 166, 215)\"],\"size\":10},\"mode\":\"markers\",\"name\":\"NCC_DMI-HIRHAM5\",\"x\":[0.4794859435892186],\"y\":[-3.769711923097435],\"type\":\"scatter\"},{\"customdata\":[[100,\"NCC_REMO2015\"]],\"hovertemplate\":\"\\u003cb\\u003eModel:\\u003c\\u002fb\\u003e NCC_REMO2015\\u003cbr\\u003e\\u003cb\\u003eTemp Bias:\\u003c\\u002fb\\u003e 0.90 \\u00b0C\\u003cbr\\u003e\\u003cb\\u003ePrec Bias:\\u003c\\u002fb\\u003e -3.9 %\\u003cextra\\u003e\\u003c\\u002fextra\\u003e\",\"marker\":{\"color\":[\"rgb(146, 150, 49)\"],\"size\":10},\"mode\":\"markers\",\"name\":\"NCC_REMO2015\",\"x\":[0.9010910335979359],\"y\":[-3.9392930754865247],\"type\":\"scatter\"},{\"customdata\":[[100,\"CNRM-CERFACS_KNMI\"]],\"hovertemplate\":\"\\u003cb\\u003eModel:\\u003c\\u002fb\\u003e CNRM-CERFACS_KNMI\\u003cbr\\u003e\\u003cb\\u003eTemp Bias:\\u003c\\u002fb\\u003e -2.6 \\u00b0C\\u003cbr\\u003e\\u003cb\\u003ePrec Bias:\\u003c\\u002fb\\u003e 50.8 %\\u003cextra\\u003e\\u003c\\u002fextra\\u003e\",\"marker\":{\"color\":[\"rgb(200, 246, 218)\"],\"size\":10},\"mode\":\"markers\",\"name\":\"CNRM-CERFACS_KNMI\",\"x\":[-2.621824679647758],\"y\":[50.86273542936121],\"type\":\"scatter\"},{\"customdata\":[[100,\"CNRM-CERFACS_ALADIN63\"]],\"hovertemplate\":\"\\u003cb\\u003eModel:\\u003c\\u002fb\\u003e CNRM-CERFACS_ALADIN63\\u003cbr\\u003e\\u003cb\\u003eTemp Bias:\\u003c\\u002fb\\u003e -1.7 \\u00b0C\\u003cbr\\u003e\\u003cb\\u003ePrec Bias:\\u003c\\u002fb\\u003e 70.3 %\\u003cextra\\u003e\\u003c\\u002fextra\\u003e\",\"marker\":{\"color\":[\"rgb(233, 25, 148)\"],\"size\":10},\"mode\":\"markers\",\"name\":\"CNRM-CERFACS_ALADIN63\",\"x\":[-1.7337617843174475],\"y\":[70.39465715222727],\"type\":\"scatter\"},{\"customdata\":[[100,\"IPSL_SMHI-RCA4\"]],\"hovertemplate\":\"\\u003cb\\u003eModel:\\u003c\\u002fb\\u003e IPSL_SMHI-RCA4\\u003cbr\\u003e\\u003cb\\u003eTemp Bias:\\u003c\\u002fb\\u003e -0.4 \\u00b0C\\u003cbr\\u003e\\u003cb\\u003ePrec Bias:\\u003c\\u002fb\\u003e 6.12 %\\u003cextra\\u003e\\u003c\\u002fextra\\u003e\",\"marker\":{\"color\":[\"rgb(82, 251, 165)\"],\"size\":10},\"mode\":\"markers\",\"name\":\"IPSL_SMHI-RCA4\",\"x\":[-0.46861886151068016],\"y\":[6.128590937718343],\"type\":\"scatter\"},{\"customdata\":[[100,\"IPSL_WRF381P\"]],\"hovertemplate\":\"\\u003cb\\u003eModel:\\u003c\\u002fb\\u003e IPSL_WRF381P\\u003cbr\\u003e\\u003cb\\u003eTemp Bias:\\u003c\\u002fb\\u003e -0.8 \\u00b0C\\u003cbr\\u003e\\u003cb\\u003ePrec Bias:\\u003c\\u002fb\\u003e 32.7 %\\u003cextra\\u003e\\u003c\\u002fextra\\u003e\",\"marker\":{\"color\":[\"rgb(190, 225, 95)\"],\"size\":10},\"mode\":\"markers\",\"name\":\"IPSL_WRF381P\",\"x\":[-0.8971256678551027],\"y\":[32.75942321691502],\"type\":\"scatter\"},{\"customdata\":[[100,\"MPI_CLMcom\"]],\"hovertemplate\":\"\\u003cb\\u003eModel:\\u003c\\u002fb\\u003e MPI_CLMcom\\u003cbr\\u003e\\u003cb\\u003eTemp Bias:\\u003c\\u002fb\\u003e -0.0 \\u00b0C\\u003cbr\\u003e\\u003cb\\u003ePrec Bias:\\u003c\\u002fb\\u003e 54.3 %\\u003cextra\\u003e\\u003c\\u002fextra\\u003e\",\"marker\":{\"color\":[\"rgb(164, 155, 82)\"],\"size\":10},\"mode\":\"markers\",\"name\":\"MPI_CLMcom\",\"x\":[-0.012765430399688908],\"y\":[54.34185803902375],\"type\":\"scatter\"},{\"customdata\":[[100,\"MPI_SMHI-RCA4\"]],\"hovertemplate\":\"\\u003cb\\u003eModel:\\u003c\\u002fb\\u003e MPI_SMHI-RCA4\\u003cbr\\u003e\\u003cb\\u003eTemp Bias:\\u003c\\u002fb\\u003e 0.04 \\u00b0C\\u003cbr\\u003e\\u003cb\\u003ePrec Bias:\\u003c\\u002fb\\u003e 16.4 %\\u003cextra\\u003e\\u003c\\u002fextra\\u003e\",\"marker\":{\"color\":[\"rgb(145, 121, 189)\"],\"size\":10},\"mode\":\"markers\",\"name\":\"MPI_SMHI-RCA4\",\"x\":[0.04996812031626846],\"y\":[16.45607454101971],\"type\":\"scatter\"},{\"customdata\":[[100,\"MPI_REMO2009\"]],\"hovertemplate\":\"\\u003cb\\u003eModel:\\u003c\\u002fb\\u003e MPI_REMO2009\\u003cbr\\u003e\\u003cb\\u003eTemp Bias:\\u003c\\u002fb\\u003e 0.85 \\u00b0C\\u003cbr\\u003e\\u003cb\\u003ePrec Bias:\\u003c\\u002fb\\u003e 4.23 %\\u003cextra\\u003e\\u003c\\u002fextra\\u003e\",\"marker\":{\"color\":[\"rgb(233, 108, 103)\"],\"size\":10},\"mode\":\"markers\",\"name\":\"MPI_REMO2009\",\"x\":[0.8586443612603661],\"y\":[4.236100432615828],\"type\":\"scatter\"},{\"hoverinfo\":\"skip\",\"marker\":{\"color\":\"black\",\"size\":10,\"symbol\":\"square\"},\"mode\":\"markers\",\"name\":\"Observations\",\"x\":[0],\"y\":[0],\"type\":\"scatter\"},{\"line\":{\"color\":\"blue\",\"width\":1.5},\"mode\":\"lines\",\"name\":\"0 Prec bias\",\"showlegend\":true,\"x\":[-2.756673228673225,1.5533013323294933],\"y\":[0,0],\"type\":\"scatter\"},{\"line\":{\"color\":\"red\",\"width\":1.5},\"mode\":\"lines\",\"name\":\"0 Temp bias\",\"showlegend\":true,\"x\":[0,0],\"y\":[-14.556739887662903,70.39465715222727],\"type\":\"scatter\"}],                        {\"template\":{\"data\":{\"histogram2dcontour\":[{\"type\":\"histogram2dcontour\",\"colorbar\":{\"outlinewidth\":0,\"ticks\":\"\"},\"colorscale\":[[0.0,\"#0d0887\"],[0.1111111111111111,\"#46039f\"],[0.2222222222222222,\"#7201a8\"],[0.3333333333333333,\"#9c179e\"],[0.4444444444444444,\"#bd3786\"],[0.5555555555555556,\"#d8576b\"],[0.6666666666666666,\"#ed7953\"],[0.7777777777777778,\"#fb9f3a\"],[0.8888888888888888,\"#fdca26\"],[1.0,\"#f0f921\"]]}],\"choropleth\":[{\"type\":\"choropleth\",\"colorbar\":{\"outlinewidth\":0,\"ticks\":\"\"}}],\"histogram2d\":[{\"type\":\"histogram2d\",\"colorbar\":{\"outlinewidth\":0,\"ticks\":\"\"},\"colorscale\":[[0.0,\"#0d0887\"],[0.1111111111111111,\"#46039f\"],[0.2222222222222222,\"#7201a8\"],[0.3333333333333333,\"#9c179e\"],[0.4444444444444444,\"#bd3786\"],[0.5555555555555556,\"#d8576b\"],[0.6666666666666666,\"#ed7953\"],[0.7777777777777778,\"#fb9f3a\"],[0.8888888888888888,\"#fdca26\"],[1.0,\"#f0f921\"]]}],\"heatmap\":[{\"type\":\"heatmap\",\"colorbar\":{\"outlinewidth\":0,\"ticks\":\"\"},\"colorscale\":[[0.0,\"#0d0887\"],[0.1111111111111111,\"#46039f\"],[0.2222222222222222,\"#7201a8\"],[0.3333333333333333,\"#9c179e\"],[0.4444444444444444,\"#bd3786\"],[0.5555555555555556,\"#d8576b\"],[0.6666666666666666,\"#ed7953\"],[0.7777777777777778,\"#fb9f3a\"],[0.8888888888888888,\"#fdca26\"],[1.0,\"#f0f921\"]]}],\"heatmapgl\":[{\"type\":\"heatmapgl\",\"colorbar\":{\"outlinewidth\":0,\"ticks\":\"\"},\"colorscale\":[[0.0,\"#0d0887\"],[0.1111111111111111,\"#46039f\"],[0.2222222222222222,\"#7201a8\"],[0.3333333333333333,\"#9c179e\"],[0.4444444444444444,\"#bd3786\"],[0.5555555555555556,\"#d8576b\"],[0.6666666666666666,\"#ed7953\"],[0.7777777777777778,\"#fb9f3a\"],[0.8888888888888888,\"#fdca26\"],[1.0,\"#f0f921\"]]}],\"contourcarpet\":[{\"type\":\"contourcarpet\",\"colorbar\":{\"outlinewidth\":0,\"ticks\":\"\"}}],\"contour\":[{\"type\":\"contour\",\"colorbar\":{\"outlinewidth\":0,\"ticks\":\"\"},\"colorscale\":[[0.0,\"#0d0887\"],[0.1111111111111111,\"#46039f\"],[0.2222222222222222,\"#7201a8\"],[0.3333333333333333,\"#9c179e\"],[0.4444444444444444,\"#bd3786\"],[0.5555555555555556,\"#d8576b\"],[0.6666666666666666,\"#ed7953\"],[0.7777777777777778,\"#fb9f3a\"],[0.8888888888888888,\"#fdca26\"],[1.0,\"#f0f921\"]]}],\"surface\":[{\"type\":\"surface\",\"colorbar\":{\"outlinewidth\":0,\"ticks\":\"\"},\"colorscale\":[[0.0,\"#0d0887\"],[0.1111111111111111,\"#46039f\"],[0.2222222222222222,\"#7201a8\"],[0.3333333333333333,\"#9c179e\"],[0.4444444444444444,\"#bd3786\"],[0.5555555555555556,\"#d8576b\"],[0.6666666666666666,\"#ed7953\"],[0.7777777777777778,\"#fb9f3a\"],[0.8888888888888888,\"#fdca26\"],[1.0,\"#f0f921\"]]}],\"mesh3d\":[{\"type\":\"mesh3d\",\"colorbar\":{\"outlinewidth\":0,\"ticks\":\"\"}}],\"scatter\":[{\"fillpattern\":{\"fillmode\":\"overlay\",\"size\":10,\"solidity\":0.2},\"type\":\"scatter\"}],\"parcoords\":[{\"type\":\"parcoords\",\"line\":{\"colorbar\":{\"outlinewidth\":0,\"ticks\":\"\"}}}],\"scatterpolargl\":[{\"type\":\"scatterpolargl\",\"marker\":{\"colorbar\":{\"outlinewidth\":0,\"ticks\":\"\"}}}],\"bar\":[{\"error_x\":{\"color\":\"#2a3f5f\"},\"error_y\":{\"color\":\"#2a3f5f\"},\"marker\":{\"line\":{\"color\":\"#E5ECF6\",\"width\":0.5},\"pattern\":{\"fillmode\":\"overlay\",\"size\":10,\"solidity\":0.2}},\"type\":\"bar\"}],\"scattergeo\":[{\"type\":\"scattergeo\",\"marker\":{\"colorbar\":{\"outlinewidth\":0,\"ticks\":\"\"}}}],\"scatterpolar\":[{\"type\":\"scatterpolar\",\"marker\":{\"colorbar\":{\"outlinewidth\":0,\"ticks\":\"\"}}}],\"histogram\":[{\"marker\":{\"pattern\":{\"fillmode\":\"overlay\",\"size\":10,\"solidity\":0.2}},\"type\":\"histogram\"}],\"scattergl\":[{\"type\":\"scattergl\",\"marker\":{\"colorbar\":{\"outlinewidth\":0,\"ticks\":\"\"}}}],\"scatter3d\":[{\"type\":\"scatter3d\",\"line\":{\"colorbar\":{\"outlinewidth\":0,\"ticks\":\"\"}},\"marker\":{\"colorbar\":{\"outlinewidth\":0,\"ticks\":\"\"}}}],\"scattermapbox\":[{\"type\":\"scattermapbox\",\"marker\":{\"colorbar\":{\"outlinewidth\":0,\"ticks\":\"\"}}}],\"scatterternary\":[{\"type\":\"scatterternary\",\"marker\":{\"colorbar\":{\"outlinewidth\":0,\"ticks\":\"\"}}}],\"scattercarpet\":[{\"type\":\"scattercarpet\",\"marker\":{\"colorbar\":{\"outlinewidth\":0,\"ticks\":\"\"}}}],\"carpet\":[{\"aaxis\":{\"endlinecolor\":\"#2a3f5f\",\"gridcolor\":\"white\",\"linecolor\":\"white\",\"minorgridcolor\":\"white\",\"startlinecolor\":\"#2a3f5f\"},\"baxis\":{\"endlinecolor\":\"#2a3f5f\",\"gridcolor\":\"white\",\"linecolor\":\"white\",\"minorgridcolor\":\"white\",\"startlinecolor\":\"#2a3f5f\"},\"type\":\"carpet\"}],\"table\":[{\"cells\":{\"fill\":{\"color\":\"#EBF0F8\"},\"line\":{\"color\":\"white\"}},\"header\":{\"fill\":{\"color\":\"#C8D4E3\"},\"line\":{\"color\":\"white\"}},\"type\":\"table\"}],\"barpolar\":[{\"marker\":{\"line\":{\"color\":\"#E5ECF6\",\"width\":0.5},\"pattern\":{\"fillmode\":\"overlay\",\"size\":10,\"solidity\":0.2}},\"type\":\"barpolar\"}],\"pie\":[{\"automargin\":true,\"type\":\"pie\"}]},\"layout\":{\"autotypenumbers\":\"strict\",\"colorway\":[\"#636efa\",\"#EF553B\",\"#00cc96\",\"#ab63fa\",\"#FFA15A\",\"#19d3f3\",\"#FF6692\",\"#B6E880\",\"#FF97FF\",\"#FECB52\"],\"font\":{\"color\":\"#2a3f5f\"},\"hovermode\":\"closest\",\"hoverlabel\":{\"align\":\"left\"},\"paper_bgcolor\":\"white\",\"plot_bgcolor\":\"#E5ECF6\",\"polar\":{\"bgcolor\":\"#E5ECF6\",\"angularaxis\":{\"gridcolor\":\"white\",\"linecolor\":\"white\",\"ticks\":\"\"},\"radialaxis\":{\"gridcolor\":\"white\",\"linecolor\":\"white\",\"ticks\":\"\"}},\"ternary\":{\"bgcolor\":\"#E5ECF6\",\"aaxis\":{\"gridcolor\":\"white\",\"linecolor\":\"white\",\"ticks\":\"\"},\"baxis\":{\"gridcolor\":\"white\",\"linecolor\":\"white\",\"ticks\":\"\"},\"caxis\":{\"gridcolor\":\"white\",\"linecolor\":\"white\",\"ticks\":\"\"}},\"coloraxis\":{\"colorbar\":{\"outlinewidth\":0,\"ticks\":\"\"}},\"colorscale\":{\"sequential\":[[0.0,\"#0d0887\"],[0.1111111111111111,\"#46039f\"],[0.2222222222222222,\"#7201a8\"],[0.3333333333333333,\"#9c179e\"],[0.4444444444444444,\"#bd3786\"],[0.5555555555555556,\"#d8576b\"],[0.6666666666666666,\"#ed7953\"],[0.7777777777777778,\"#fb9f3a\"],[0.8888888888888888,\"#fdca26\"],[1.0,\"#f0f921\"]],\"sequentialminus\":[[0.0,\"#0d0887\"],[0.1111111111111111,\"#46039f\"],[0.2222222222222222,\"#7201a8\"],[0.3333333333333333,\"#9c179e\"],[0.4444444444444444,\"#bd3786\"],[0.5555555555555556,\"#d8576b\"],[0.6666666666666666,\"#ed7953\"],[0.7777777777777778,\"#fb9f3a\"],[0.8888888888888888,\"#fdca26\"],[1.0,\"#f0f921\"]],\"diverging\":[[0,\"#8e0152\"],[0.1,\"#c51b7d\"],[0.2,\"#de77ae\"],[0.3,\"#f1b6da\"],[0.4,\"#fde0ef\"],[0.5,\"#f7f7f7\"],[0.6,\"#e6f5d0\"],[0.7,\"#b8e186\"],[0.8,\"#7fbc41\"],[0.9,\"#4d9221\"],[1,\"#276419\"]]},\"xaxis\":{\"gridcolor\":\"white\",\"linecolor\":\"white\",\"ticks\":\"\",\"title\":{\"standoff\":15},\"zerolinecolor\":\"white\",\"automargin\":true,\"zerolinewidth\":2},\"yaxis\":{\"gridcolor\":\"white\",\"linecolor\":\"white\",\"ticks\":\"\",\"title\":{\"standoff\":15},\"zerolinecolor\":\"white\",\"automargin\":true,\"zerolinewidth\":2},\"scene\":{\"xaxis\":{\"backgroundcolor\":\"#E5ECF6\",\"gridcolor\":\"white\",\"linecolor\":\"white\",\"showbackground\":true,\"ticks\":\"\",\"zerolinecolor\":\"white\",\"gridwidth\":2},\"yaxis\":{\"backgroundcolor\":\"#E5ECF6\",\"gridcolor\":\"white\",\"linecolor\":\"white\",\"showbackground\":true,\"ticks\":\"\",\"zerolinecolor\":\"white\",\"gridwidth\":2},\"zaxis\":{\"backgroundcolor\":\"#E5ECF6\",\"gridcolor\":\"white\",\"linecolor\":\"white\",\"showbackground\":true,\"ticks\":\"\",\"zerolinecolor\":\"white\",\"gridwidth\":2}},\"shapedefaults\":{\"line\":{\"color\":\"#2a3f5f\"}},\"annotationdefaults\":{\"arrowcolor\":\"#2a3f5f\",\"arrowhead\":0,\"arrowwidth\":1},\"geo\":{\"bgcolor\":\"white\",\"landcolor\":\"#E5ECF6\",\"subunitcolor\":\"white\",\"showland\":true,\"showlakes\":true,\"lakecolor\":\"white\"},\"title\":{\"x\":0.05},\"mapbox\":{\"style\":\"light\"}}},\"legend\":{\"font\":{\"size\":8},\"x\":1.02,\"y\":0.95},\"margin\":{\"l\":40,\"r\":40,\"t\":40,\"b\":40},\"title\":{\"text\":\"Friuli-Venezia Giulia\"},\"xaxis\":{\"title\":{\"text\":\"Mean temperature bias (\\u00b0C)\"}},\"yaxis\":{\"title\":{\"text\":\"Mean precipitation bias (%)\"}},\"height\":400,\"width\":600},                        {\"responsive\": true}                    ).then(function(){\n",
       "                            \n",
       "var gd = document.getElementById('58ab1288-036b-4531-8e77-fcbfba46fbde');\n",
       "var x = new MutationObserver(function (mutations, observer) {{\n",
       "        var display = window.getComputedStyle(gd).display;\n",
       "        if (!display || display === 'none') {{\n",
       "            console.log([gd, 'removed!']);\n",
       "            Plotly.purge(gd);\n",
       "            observer.disconnect();\n",
       "        }}\n",
       "}});\n",
       "\n",
       "// Listen for the removal of the full notebook cells\n",
       "var notebookContainer = gd.closest('#notebook-container');\n",
       "if (notebookContainer) {{\n",
       "    x.observe(notebookContainer, {childList: true});\n",
       "}}\n",
       "\n",
       "// Listen for the clearing of the current output cell\n",
       "var outputEl = gd.closest('.output');\n",
       "if (outputEl) {{\n",
       "    x.observe(outputEl, {childList: true});\n",
       "}}\n",
       "\n",
       "                        })                };                });            </script>        </div>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "plots_list_plotly[100]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 16,
   "id": "3d5517c0-10c8-4367-ad44-bee9fef53054",
   "metadata": {
    "tags": []
   },
   "outputs": [
    {
     "data": {
      "text/html": [
       "\n",
       "        <iframe\n",
       "            width=\"100%\"\n",
       "            height=\"650\"\n",
       "            src=\"http://127.0.0.1:8050/\"\n",
       "            frameborder=\"0\"\n",
       "            allowfullscreen\n",
       "            \n",
       "        ></iframe>\n",
       "        "
      ],
      "text/plain": [
       "<IPython.lib.display.IFrame at 0x2699e1bd250>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "# Initialize Dash app\n",
    "app = Dash(__name__)\n",
    "\n",
    "# App layout (select region and show main graph + details)\n",
    "app.layout = html.Div([\n",
    "    html.H4(\"EURO-CORDEX Climate Model Biases and Uncertainty at NUTS2 level\"),\n",
    "    dcc.Dropdown(\n",
    "        id='region-select',\n",
    "        options=[{'label': name, 'value': i} for i, name in enumerate(regions_ex2['NUTS_NAME'])],\n",
    "        value=0\n",
    "    ),\n",
    "    dcc.Graph(id='main-graph'),\n",
    "    html.Div(id='detail-output')\n",
    "])\n",
    "\n",
    "# Main graph update on region selection\n",
    "@app.callback(\n",
    "    Output('main-graph', 'figure'),\n",
    "    Input('region-select', 'value')\n",
    ")\n",
    "def update_main_plot(region):\n",
    "    region_name = regions_ex2['NUTS_NAME'][region]\n",
    "    x_data, y_data, model_labels, colors = [], [], [], []\n",
    "\n",
    "    for model in range(combinations):\n",
    "        index = region + (len(regions_ex2) * model)\n",
    "        x_data.append(median_bias_tas[index])\n",
    "        y_data.append(median_bias_pr[index])\n",
    "        model_labels.append(label_list[model])\n",
    "        colors.append(color_list[model])\n",
    "\n",
    "    fig = go.Figure()\n",
    "\n",
    "    for i in range(combinations):\n",
    "        fig.add_trace(go.Scatter(\n",
    "            x=[x_data[i]],\n",
    "            y=[y_data[i]],\n",
    "            mode='markers',\n",
    "            marker=dict(color=colors[i], size=10),\n",
    "            name=model_labels[i],\n",
    "            customdata=[[region, model_labels[i]]],\n",
    "            hovertemplate=(\n",
    "                f\"<b>Model:</b> {model_labels[i]}<br>\"\n",
    "                f\"<b>Temp Bias:</b> {str(x_data[i])[:4]} °C<br>\"\n",
    "                f\"<b>Prec Bias:</b> {str(y_data[i])[:4]} %<extra></extra>\"\n",
    "            )\n",
    "        ))\n",
    "\n",
    "    fig.add_trace(go.Scatter(\n",
    "        x=[0], y=[0],\n",
    "        mode='markers',\n",
    "        marker=dict(color='black', size=10, symbol='square'),\n",
    "        name='Observations',\n",
    "        hoverinfo='skip'\n",
    "    ))\n",
    "\n",
    "    fig.add_trace(go.Scatter(x=[min(x_data + [0]), max(x_data + [0])], y=[0, 0],\n",
    "                             mode='lines',\n",
    "                             line=dict(color='blue', width=1.5),\n",
    "                             name='0 Prec bias'))\n",
    "\n",
    "    fig.add_trace(go.Scatter(x=[0, 0], y=[min(y_data + [0]), max(y_data + [0])],\n",
    "                             mode='lines',\n",
    "                             line=dict(color='red', width=1.5),\n",
    "                             name='0 Temp bias'))\n",
    "\n",
    "    fig.update_layout(\n",
    "        title=f\"{region_name}\",\n",
    "        xaxis_title=\"Mean temperature bias (°C)\",\n",
    "        yaxis_title=\"Mean precipitation bias (%)\",\n",
    "        legend=dict(font=dict(size=8)),\n",
    "        height=400,\n",
    "        width=600\n",
    "    )\n",
    "\n",
    "    return fig\n",
    "\n",
    "# Uncertainty plot when clicking a model point\n",
    "@app.callback(\n",
    "    Output('detail-output', 'children'),\n",
    "    Input('main-graph', 'clickData')\n",
    ")\n",
    "def display_details(clickData):\n",
    "    if clickData is None:\n",
    "        return html.Div(\"Click a model to see more details.\")\n",
    "\n",
    "    model = clickData['points'][0]['customdata'][1]\n",
    "\n",
    "    detailed_fig = go.Figure()\n",
    "    detailed_fig.add_trace(go.Scatter(\n",
    "        x=np.arange(10),\n",
    "        y=np.random.rand(10),\n",
    "        mode='lines+markers',\n",
    "        name=f'{model} Time Series'\n",
    "    ))\n",
    "    detailed_fig.update_layout(title=f\"Uncertainty for {model}\")\n",
    "\n",
    "    return dcc.Graph(figure=detailed_fig)\n",
    "\n",
    "# Run app\n",
    "if __name__ == '__main__':\n",
    "    app.run_server(debug=True)\n",
    "    \n",
    "fig.write_html(\"bias_plotly.html\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "277cc537-aed4-4810-a0d9-45a1097e631c",
   "metadata": {
    "tags": []
   },
   "outputs": [],
   "source": [
    "\n",
    "\n",
    "\n",
    "\n",
    "\n",
    "\n",
    "\n",
    "\n",
    "\n",
    "\n",
    "\n",
    "\n",
    "\n",
    "\n",
    "\n",
    "\n",
    "\n",
    "\n",
    "\n",
    "\n",
    "\n",
    "\n",
    "\n"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "paper_grano",
   "language": "python",
   "name": "paper_grano"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.11.5"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
