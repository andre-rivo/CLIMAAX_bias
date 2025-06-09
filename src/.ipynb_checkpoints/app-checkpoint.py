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
      "image/png": "iVBORw0KGgoAAAANSUhEUgAAAvQAAAGQCAYAAADBUfTHAAAAAXNSR0IArs4c6QAAIABJREFUeF7svQd8VMX6///Z3XQgEBI6SO+EjiiggEgTKaJ0sNG8FkRQVL4WVC4oNhDkKihSpSlIExAQAWleUOlFCL0lBAiE9GT/vxnu5r/Z7GY3u3tmz4bPeb3u6wI7Z56Z93MS32f2OXMMZrPZDB4kQAIkQAIkQAIkQAIkQAJ+ScBAoffLvHHQJEACJEACJEACJEACJCAJUOh5IZAACZAACZAACZAACZCAHxOg0Ptx8jh0EiABEiABEiABEiABEqDQ8xogARIgARIgARIgARIgAT8mQKH34+Rx6CRAAiRAAiRAAiRAAiRAoec1QAIkQAIkQAIkQAIkQAJ+TIBCn0fyDh07jSGjJ2HIgEcxuN8jfpxmDt1C4NuFP+ObBavxzadjULdmJYIhARIgARIgARIgAb8noGuhF/L12ddLHEKeNuFltG3RKM8kjJ04E5u3/+WWwNkTesuYvBXbXwQzLv4GBr00AWVLRULMPSw0xKsXv6X/cxdjc/RboWxJzJs6FiUii2X/uyc5tcfbk/68CoGdkQAJkAAJkAAJkIAbBHQv9J6upnoia54KvbMV/qTkFLw4dopMmxaS7Mb14PAULYVe5GjF+u0YNbx3jm9CLHx2/3Ukx2ee5JRC782rgn2RAAmQAAmQAAnogUCBF3pPIDsTcmd9OxN2T/t3Ft8fPndFzjfv+AsxZy55pezJX74R8YfccYwkQAIkQAIkQAL6IOD3Qm9ZOW4cXR0DeraXNe83E5PQvFFtueo9fvI8/Hngn+yyDUv7Xl3b2F0NLl2yOCa8OVRmx9EKfX6+NcirREd8tnTVbzlKSuyVntiW9wgJFnP6bNzzGDVuOixlKpY525bDCCG2fBMg5mVbxmL7ufWlaWkr/k2U3AjOFj6OymS6d2yZ3Savy9zCt23LRi61t/Rlmb+lFMfSz4SxQ3OUYNnLtaMVeutrxPqbAevxO+Krjx9ljoIESIAESIAESOBuJVBghF5IrT3hspU/1UJvfcNhEWFxsdn7d4tYW5ee2LupsJSoWIu5ozj2VsBdWaW2HYu9/sW/vT7+a3z01vDsGndH47D3A5af5xGsz1ch9G9MmIHhA7tlPzjr7NuWu/UXCOdNAiRAAiRAAiTgewK6F3pHD8VapNeZQPpa6EWKbccg/s0izJbVd4swWn9D4GhF2l5/oq3tir+zlWvr1XbrS9FyXu3qFbNr+51xtj7f3jcPjoQ+P992qFqhd/RjKXI2dsJMtx6w9v2POkdAAiRAAiRAAiRQUAnoXuidCZ8z0dRa6O2VqzgqabEunXF1ldki/9YimZfQW/PKS6xFH5djr+V6GNfCU8S13l0mL86Wbwysf0js7U5j+0Pk6JsCe7sbWX9r4So7d0tuLOO0N47wwmEU+oL625DzIgESIAESIAE/JUCh/1/i7K2Qe6OGXnRvK8P2RDOvOnbRh7VIuir09kTb+jp1VKJkb5tPe0Jv4VM0vHAO+c/PCr34BiavLUDtlSFpLfSWuSbcTMwh71yh99Pfchw2CZAACZAACRRwAhR6BUIvQljXssfGX89VuuGoPMbe9eeq0Lsq1rYr0rbbR9q7KbHMyfphUut+bB/2tTcPVx6K9YXQO/rmgEJfwH8bcnokQAIkQAIk4KcE7lqht60f13KFXlwbFnnt26Md9h06IS8X673nnZUOWV9frgq9bZ1+Xteopa2jHWpsx5fXQ6L5uZGwfIvgaJU+P0Jv+0Zfd0tuHPGl0PvpbzkOmwRIgARIgAQKOIG7TujtrSxbb1NoLbTeKrkRMW23QrS3Cu5Iqm13k3FV6C1zFS9tshVmEWvDlj1yu0hXXhpl74bD3kq2pe7clRp6Wy72mFj6y6uG3t7NhfWWmtbnurJtpb0bIcu/sYa+gP9G5PRIgARIgARIwA8J3JVCbyvXQtKmfDACX81dAW/vQ299TVjkNC8pdLS3e15Ca4mRV6mI9T70or31GOw9/Gnp09KuZFSxXPvQiza254pxisOVkhtrNpabJ/EOAevDla1IRXtbbuKGwrJPv/U7B1wRetGf7TMN4kavfeum3OXGD3/JccgkQAIkQAIkUNAJ6FroCzp8zo8ESIAESIAESIAESIAEPCVAofeUIM8nARIgARIgARIgARIgAR8SoND7ED5DkwAJkAAJkAAJkAAJkICnBCj0nhLk+SRAAiRAAiRAAiRAAiTgQwIUeh/CZ2gSIAESIAESIAESIAES8JQAhd5TgjyfBEiABEiABEiABEiABHxIgELvQ/gMTQIkQAIkQAIkQAIkQAKeEqDQe0qQ55MACZAACZAACZAACZCADwlQ6H0In6FJgARIgARIgARIgARIwFMCFHpPCfJ8EiABEiABEiABEiABEvAhAQq9D+EzNAmQAAmQAAmQAAmQAAl4SoBC7ylBnk8CJEACJEACJEACJEACPiRAofchfIYmARIgARIgARIgARIgAU8JUOg9JcjzSYAESIAESIAESIAESMCHBCj0PoTP0CRAAiRAAiRAAiRAAiTgKQEKvacEeT4JkAAJkAAJkAAJkAAJ+JAAhd6H8BmaBEiABEiABEiABEiABDwlQKH3lCDPJwESIAESIAESIAESIAEfEqDQ+xA+Q5MACZAACZAACZAACZCApwQo9J4S5PkkQAIkQAIkQAIkQAIk4EMCd7XQX4xPtou+bGQoHH3mw1wx9P8IhASZEBZswrVbaWSiUwLFiwQhKTUTKWmZOh3h3T2s8LBAFA4NwM2kdCQmZ9zdMHQ6+7CQAASZDLhxO12nIyw4wxL/zedBAv5OgEJvJ4MUen1f1hR6fedHjI5Cr+8cUej1nR8xOgq9uhxR6NWxZiTtCFDoKfTaXV0a9Uyh1wisF7ul0HsRpgZdUeg1gOrlLin0XgaaR3cUenWsGUk7AhR6Cr12V5dGPVPoNQLrxW4p9F6EqUFXFHoNoHq5Swq9l4FS6NUBZSSfEKDQU+h9cuF5EpRC7wk9NedS6NVwdjcKhd5dcurOo9CrY80VenWsGUk7AhR6Cr12V5dGPVPoNQLrxW4p9F6EqUFXFHoNoHq5Swq9l4FyhV4dUEbyCQEKPYXeJxeeJ0Ep9J7QU3MuhV4NZ3ejUOjdJafuPAq9OtZcoVfHmpG0I0Chp9Brd3Vp1DOFXiOwXuyWQu9FmBp0RaHXAKqXu6TQexkoV+jVAWUknxCg0FPofXLheRKUQu8JPTXnUujVcHY3CoXeXXLqzqPQq2PNFXp1rBlJOwIUegq9dleXRj1T6DUC68VuKfRehKlBVxR6DaB6uUsKvZeBcoVeHVBG8gkBCj2F3icXnidBKfSe0FNzLoVeDWd3o1Do3SWn7jwKvTrWXKFXx5qRtCNAoafQa3d1adQzhV4jsF7slkLvRZgadEWh1wCql7uk0HsZKFfo1QFlJJ8QoNBT6H1y4XkSlELvCT0151Lo1XB2NwqF3l1y6s6j0KtjzRV6dawZSTsCFHoKvXZXl0Y9U+g1AuvFbin0XoSpQVcUeg2gerlLCr2XgXKFXh1QRvIJAQo9hd4nF54nQSn0ntBTcy6FXg1nd6NQ6N0lp+48Cr061lyhV8eakbQjQKGn0Gt3dWnUM4VeI7Be7JZC70WYGnRFodcAqpe7pNB7GShX6NUBZSSfEKDQU+h9cuF5ElSp0N9OB87fgiEjC+aK4UDhIE+GftecS6HXd6op9PrOjxgdhV5djrhCr441I2lHgEJPodfu6tKoZxVCb7ieAsPne2BcfjzHLLI6V0HWa/cCJcI0ml3B6JZCr+88Uuj1nR8Kvdr8UOjV8mY0bQhQ6Cn02lxZGvaqtdAbkjNg7LsShpgbdmdhLlcYmYu7A0WDNZylf3dNodd3/ij0+s4PhV5tfij0ankzmjYEKPQUem2uLA171VrojR//AePcg3nOIKtnDWS910rDWfp31xR6feePQq/v/FDo1eaHQq+WN6NpQ4BCT6HX5srSsFethd7UfxUMB+LynIG5UlFkrnpcw1n6d9cUen3nj0Kv7/xQ6NXmh0KvljejaUOAQk+h1+bK0rBXrYU+4P75QGJa3jMwGpDx59OAyaDhTP23awq9vnMXnnwThU8cwc3a9ZEYEKrvwd6lo+NDseoST6FXx5qRtCNAoafQa3d1adSz1kJv6vojDKcT8h59ZCgyfuun0Qz9v1sKvb5zGLFlA0If747bW3cgoU5DfQ/2Lh0dhV5d4in06lgzknYEKPQUeu2uLo161lzoX/8Nhp9j8hy9+YEKyJzeXqMZ+n+3FHp955BCr+/8iNFR6NXliEKvjjUjaUeAQk+h1+7q0qhnrYVe7Dsf0GclcDPV7gzMoQHIWtgV5qoRGs3Q/7ul0Os7hxR6feeHQq82PxR6tbwZTRsCFHoKvTZXloa9ai70AAwbTsM06le7s8j64AFk9aiu4Qz9v2sKvb5zSKHXd34o9Grzo1Toj+xH1hcfAKeOAwnXgTLlYWjVHoZhrwKFi6idOKMVKAIUegq9313QKoReQknNhPHXM8DvF2DIzEJWq3JAu0oQK/Q88iZAodf3FUKh13d+KPRq86NE6BNvwfzlBJiXfgeYzbknGFUShtHjYejQXe3kGa3AEKDQU+j97mJWJvR+R0Y/A6bQ6ycX9kZCodd3fij0avOjROj/OYysvm3znti9D8L4n6VuTz4pORXjJ8/Fpdh4GGBAZPFwjB0xEAk3b+PDad/jcuw1hAQHom+PdujRKed7VD75ajEOHImByWhEWGgIRj/XG5XvKYMNW/dg5frtmPDmUBQpHIbjMecx4Yv5MJvNiIwIx5gX+snPd+09jBs3E+W5QYEBGNK/i2zbu1tbhIUGyz8fP3kOj7a/3+78Tp+7nGuMpaLulLXe37Su/H8xv/c/n4PObZuj9f0NEH/9JsaM/wpD+nXJbuM2vAJwIoWeQu93lzGFXv8po9DrO0cUen3nh0KvNj8qhN783RcwT/u304kZfzsOFCnqtJ29BknJKZi5YA2GDugixVoctxKTMH7KPIwY/DjKlY5CZmYWYs5eRPXK5XN08cW3P6Jfj3YoEVks+99T09IxZeYPyMjMQqe2zdA4ugYOHTuNg0dj0Kf7Q9i2+wD2Hz6BF555TJ4jRL97x1aoW7MSbMdifZ7t2B2N8fylO++Dadui0f+EPgUTp34PgwHyRmXLzn2Yu3Q9hgzoItucvXAFR/45i/JlSqBWtXtw7cZN7D8cg4hihdGgTjVcvZaAy3HX5L8XCy+McqVLyJuYCuVKokaVnDzcSoCPT6LQU+h9fAnmPzyFPv/MVJ9BoVdNPH/xKPT54+WL1tzlRh11JUL/7kswr17idFLG2WuA6KZO2zkS+rcnzUJoSDACTCYUDS+EFs2isWnbHinAeR0ffD4XiUnJCA0ORnBwIJ7s1VFK+brNf6Bls2hs3bUPI4c+kS30YuV95S87kJmZiZ6PPGhX6K3Hcj3hFprUryH7tT2E7K9Y/3uuMW7e8VcuoRc3LOGFw1CqRHGcOH1efqNQp2YlBAcFYtO2PzG43yM4dPw0qlYsh9lL1uKp3p2w5++jCAgIQPGIIliy8je8MqwX5ixZh5DgIPTv+TC+/X4NRg3vjaji7t1IuZUsDU6i0FPoNbistO2SQq8tX2/0TqH3BkXt+qDQa8fWWz1T6L1F0nk/SoT+36/BvGyu08F4KvTT56xA765tpNQbjUZcuHQVqzfuyCXLi1f8ilmL1qJ+nSp479VnMWP+KnRqe68sozEYDCgaXhjL1mxBeJFCqF29ovxcSO+VuOuYPPMHmExGlIwqhndGPSVvHsRhu0JvPZZjJ8/h3MVY+S2Ap0L/WOdWGPfJbAzq1QExZy6hSsUyOHHqAu5rUhfRtSrL7q2/EYiLv4Gv5q5E8ya1ce5CnJT+bxf+LM8TK/vW43aaIB03oNBT6HV8edofGoVe/ymj0Os7RxR6fedHjI5Cry5HSoT+hzkwTxyT96QMBhh/jwFCwtyavL2SG1FmMnHqAin0QtbFIWrqxeq99WFbciPKYETd/X2N6yAoKBB79x+Twizq2kXJTc8urfHx9IVo2qAWOrS+842CuyU3jsb458Hjsl/rkhtLSdG5i3GoWL4UFizbKMX83IVYVL6nLB5oHi3POXriLHbsOYhn+z4ibyQW/fQrmjSogVNnL1Po3bq6dH7SxfhkuyMUP9yOPtP5lO6K4VHo9Z9mCr2+c0Sh13d+KPRq86NC6HHzBrKeaAXE36kLt3cY+g2F4dXxbk9eCL1tyc3g/l2w79AJTJu1HDWqVkBc/HU0a1hbSq31YVty07RBTZw6ewnDBnaVzf45dR6LV2xG1w4tcPSfM7KGXkj/uE9n4+nenRBdu4rbQi/637Z7f64xClFfuHwTypaKkmPo1rGFrNu3fkbAstJes0oF+cBsicgIpKamYdRzfTBv6XpcuXodaenpch5ipV6s6HOF3u1LTL8nUuj1m5u8Rkah13/eKPT6zhGFXt/5odCrzY8SoRdT2r0FWS/0sb9tZZ2GMH63GggI1GzyYpeY0JAgWVLjq+NmYpIUbuujUFhI9kO83h6juOkICQlGYMCdsqCCfLDkxk52uUKv70ueQq/v/IjRUej1nSMKvb7zQ6FXmx9lQi+mJfajnzcd5h2/AudOAY3ug/GxAcCDuR8WVUtBTTRnQq9mFAUzCoWeQu93VzaFXv8po9DrO0cUen3nh0KvNj9KhV7t1BjtLiJAoafQ+93lTqHXf8oo9PrOEYVe3/mh0KvND4VeLW9G04YAhZ5Cr82VpWGvFHoN4Xqpawq9l0Bq1A2FXiOwXuyWu9x4EaaTrij06lgzknYEKPQUeu2uLo16ptBrBNaL3VLovQhTg64o9BpA9XKXFHovA82jOwq9OtaMpB0BCj2FXrurS6OeKfQagfVitxR6L8LUoCsKvQZQvdwlhd7LQCn06oAykk8IUOgp9D658DwJSqH3hJ6acyn0aji7G4VC7y45dedR6NWx5gq9OtaMpB0BCj2FXrurS6OeKfQagfVitxR6L8LUoCsKvQZQvdwlhd7LQLlCrw4oI/mEAIWeQu+TC8+ToBR6T+ipOZdCr4azu1Eo9O6SU3cehV4da67Qq2PNSNoRoNBT6LW7ujTqmUKvEVgvdkuh9yJMDbqi0GsA1ctdUui9DJQr9OqAMpJPCFDoKfQ+ufA8CUqh94SemnMp9Go4uxuFQu8uOXXnUejVsVa5Qr8vIwbvJs7D8cwLuJZ1CxVMJdAhqAleL9QL4YYwdZNmpAJHgEJPofe7i5pCr/+UUej1nSMKvb7zI0ZHoVeXIxVCf9OchPG3v8c3yethhjnX5EoZi2Fi4WfxWHALdRNnpAJFgEJPofe7C5pCr/+UUej1nSMKvb7zQ6FXmx8VQn8o8wxaXRud58RaB0bjp2Lvuj35pORUjJ88F5di42GAAZHFwzF2xEAk3LyND6d9j8ux1xASHIi+PdqhR6dWOeJ88tViHDgSA5PRiLDQEIx+rjcq31MGG7buwcr12zHhzaEoUjgMx2POY8IX82E2mxEZEY4xL/STn+/aexg3bibKc4MCAzCkfxfZtne3tggLDZZ/Pn7yHB5tf3+u+Ym+lq/dhsUrNiM5NQ1lS0Xi7VeexNnzV7B09RaMf/1Z2e+yn7di595DeO/VZ7F4xa/o1rGlHEP89ZtyDD06t8J7n87BxStX8caL/dE4ukZ2LGs2GRlZeO7JbmjZrJ6cx7rNf2DuD78gOTkVUZFF8eZLA1C1Ytlcc7d0dvDYKSxZuRlvvDhAzk0vB4WeQq+Xa9HlcVDoXUbls4YUep+hdykwhd4lTD5txBV6dfhVCP3nScvx/u0FTid1OmoOihoKOW1nr0FScgpmLliDoQO6SAEWx63EJIyfMg8jBj+OcqWjkJmZhZizF1G9cvkcXXzx7Y/o16MdSkQWy/731LR0TJn5AzIys9CpbTMpyIeOncbBozHo0/0hbNt9APsPn8ALzzwmzxGi371jK9StWQm2Y7E+z3bsJ89cxNyl6/HWyCcRGGBCckoazOYs7P7rCBYu34Qne3VEgzpVMXHqAiSnpOLfbwzBN9+vyR5vXPwNLPxpk5yjozjW4xFy/+lXS/DGS/1x6uwlrNm4E2Ne6C9ji5uf5NRURBQtkmvuYtwXL1/Fl7N/QkhwkLzpsXB2K2FePolCT6H38iWlfXcUeu0ZexqBQu8pQW3Pp9Bry9cbvVPovUHRtT5UCP3zt6ZhYcpvTgf0S8QENAv4/1eWnZ5g1UBI69uTZiE0JBgBJhOKhhdCi2bR2LRtj1ypz+v44PO5SExKRmhwMIKDA6VEi/7E6nXLZtHYumsfRg59IluYxcr7yl92IDMzEz0fedCu0FuP5XrCLTSpX0P2a3vcSEjEyHenoU2LhmjVLBoVK5SWcr15x1+Ii0/AydMX5E2CuBn5799H8dbIQVLILeMVAl44LFSu6rsi9JfjruP7ZRvx2vN9Mf/HDahSsQzatmiUY1j/nDqfa+7i5mj+sg14rPMD8tsE6xun/ORJq7YFQuhF0l8cO0UyqlC2JOZNHZt9l/ntwp/x2ddL5GfNG9XGtAkvZ99RXYxPtstV/HA7+kyrRLBf1wlQ6F1n5auWFHpfkXctLoXeNU6+bEWhV0dfhdC/cutrzE7Z4HRSngr99Dkr0LtrGyn1RqMRFy5dxeqNO3IJvShZmbVoLerXqSJLWGbMX4VObe+VJSwGgwFFwwtj2ZotCC9SCLWrV5SfjxreG1firmPyzB9gMhlRMqoY3hn1lLx5EIftCr31WI6dPIdzF2Plqrq9Q3wbsO/QSezYc1CuzI9/fTDOXrgim169loCtu/Zj9PDechxC6GfMX509XlFyI248rG84xDcI1nN8/YX+coXfaDTg7PlYTP7gJZQpWRzCEe0JvTjXdu6rN+5Ew7rVULxYESxYthHP9O2M0iWKO82pqgZ+L/RC5j+atjCHxFvg2X42duLMOxfdm0Pl/1PoVV1m3o1DofcuTy16o9BrQdV7fVLovcdSq54o9FqRzd2vCqGflfwLRifOyHNSou79fIn5CIN7ddn2Sm6EDAuRFSv0QtbFIcpKxOq99WFbciNWo0Xd/X2N6yAoKBB79x/DfU3qolRUhCy56dmlNT6evhBNG9RCh9ZN7Qq9dflPXiU3osTGZDTIOOKYuWA1alW7BxmZmfLvjevVQMKtRHmT8vmMpVLoPSm5+f2Pg9iz7yhee76frP3f/ddhebMgbmTSMzJx7XoCvvh2Wa65R0WE41LsNaSlpWPjtr0Y+Hh7NGtYS92F6iSSXwu9qJt6/s3JGDf6afl1jO0hBL5qpXIY3O8R+ZGt4FPodXMd5msgFPp84fJJYwq9T7C7HJRC7zIqnzWk0KtDr0Lob5gT0fzaSMRm3XA4sedCu2Bi4Wfcnri9kpvB/btg36ETmDZrOWpUrYC4+Oto1rB2thdZgtmW3DRtUFPWlw8b2FU2ESUoosyka4cWOPrPGVlDL6R/3Kez8XTvToiuXcXtGvrT5y5jzPivcE+5UjKWEPcxz/fFnv3H5N8t5TDC+fIS+qf7dMYn/1mEU2cvY+yIAfKbBcthfbMj+p+9eJ38aMDj7TF99k/Yf/gkypUpIb9FePTh+3Htxs1ccxc19+LbCHs3Tm4nzYsn+rXQizu+IaMn4eb/u6gsR/eOLeUKvAAuynBa3hudfeGK9qPGfYnPxr0gbwAo9F68khR2RaFXCNvNUBR6N8EpOo1Crwi0B2Eo9B7Ay+epKoReDOm39P3oeeMDu9tWNgqoivUR/0YgAvI5etebi4dBQ0OC5Eq0rw7ha6mpaTnCFwoLkaXQouxGHMH/W6lXOUZRn5+WniH5+Ovh10IvVtznLf0luy5e3L0NemkCenVtg349HpJCP6hXh+y7O1uhv5l05+KxPcLDAuHoM39NdEEad4DJiKAAA5JS73wdx0N/BMKCTUjLMMvdEfR8ZGWZZU3l3XaEbViHgB7dkLpjJ1IbNLnbpu8X8w0MMCLAaEByGn/PaZ0w8d98VYfYj35a0kpsTPsLMZmX0SKwNgaFPozOQXfKVgr6kZfQF/S5az2/AiX0ApZ4wGH7Hwcw6e3hGPPB13mu0CcmZ9jlWzg0AI4+0zoh7N85gQCTAULqU/gfOuewfNRCfIsiZD4jM/cLVHw0JLth71ahD9mwFgHduyFtxy6kNaTQ6+matIxFCL2410xN1/dNsR7Z5XdM4r/5PEjA3wn4tdCLFXdRvzV94sgcu9qILY5E2Q1r6P398rQ/fpbc6D+vLLnRd45YcqPv/IjRseRGXY5UldyomxEj3Y0E/FroLXXypUsWlwJvKbl5/cV+ssyGu9wUzEuaQq//vFLo9Z0jCr2+80OhV5sfCr1a3oymDQG/FnqBxCLx4slkcYh9Ui272oi/cx96bS4cX/ZKofclfddiU+hd4+SrVhR6X5F3PS5X6F1n5WlLCr2nBHm+Hgg4FHprEbYdqK0062Ei7oyBu9y4Q83351DofZ8DZyOg0Dsj5NvPKfS+5e9KdAq9K5S804ZC7x2O7MW3BHIJvag7X7F+e643rlqGab0ibtki0rdTcD86hd59dr48k0LvS/quxabQu8bJV60o9L4i73pcCr3rrDxtSaH3lCDP1wOBbKG3iHrZUpHZ20DmNUBL/frFK/F239Kqh8k5GwOF3hkhfX5OoddnXqxHRaHXd44o9PrOjxgdhV5djij06lgzknYEcgj97r+OyDdk5edYvXEnmjeqnb3LTH7O9XVbCr2vM+BefAq9e9xUnkWhV0k7/7Eo9PlnpvoMCr064hR6dawZSTsCfv9QrCdoKPSe0PPduRR637F3NTKF3lVSvmlHofcN9/xEpdDnh5ZnbSn0nvHj2fogQKG3kwfxw+1I9vWRtrt7FBR6/eefQq/vHFHo9Z0fMToKvbocqRT6uLgM7Np+G9evZSIlJQtFiphQsVIQmjYPQ1DQ3ffWanWyoM1dAAAgAElEQVRZLviRnAq9eHnTkNGTIF7XK44KZUv6bc28bTq5Qu+fFziFXv95o9DrO0cUen3nh0KvNj8qhD4tzYzdO2/j4P4Uu5MLCzOi5YOFUK16sNrJM1qBIZCn0IsHX9+YMAPDB3ZD3ZqV5KTFy5rmLf3FpQdn9U6JQq/3DNkfH4Ve/3mj0Os7RxR6feeHQq82PyqEPj4+E0u+v57nxMpXCETXHkXdnnxScirGT56LS7HxMMCAyOLhGDtiIBJu3saH077H5dhrCAkORN8e7dCjU6sccT75ajEOHImByWhEWGgIRj/XG5XvKYMNW/dg5frt8uWdRQqH4XjMeUz4Yj7MZjMiI8Ix5oV+8vNdew/jxs1EeW5QYACG9O8i2/bu1hZhocHyz8dPnsOj7XM/pynGvWTlZtk2OCgQi1b8iqqVysKcZcaMBaswcewwlC5RHN8tWotuHVvKfpau3oLxrz8r4y37eSt27j2EkUN74ePpi3DxylW88WJ/NI6u4TZLfzwxT6EXO9+8Pv5rfPTW8OyHXu39mz9OXIyZQu+fmaPQ6z9vFHp95igrIxVJSbGI2LYZJZ56Dre37kBCnYb6HOxdPiqW3Ki7AFQI/Z97k7F7x22nk3p2WCSCg90rvRGLsDMXrMHQAV2k6IrjVmISxk+ZhxGDH0e50lHIzMxCzNmLqF65fI6xfPHtj+jXo12ODU5S09IxZeYPyMjMQqe2zaQgi6qNg0dj0Kf7Q9i2+wD2Hz6BF555TPYlRL97x1ZyAdh2LNbn2UKwtBUvBV2+dhsCAkzo2/0h/Lbzb2zdtR9V7imDQU90gGWMB4+dwsLlm/Bkr45oUKcqJk5dgOSUVPz7jSE4dfZy9vicwi5gDXIIvZD1z2csxVsjB8mLgSv0BSzbBWQ6FHr9J5JCr78cnTm2EqeO/ACzOQtl919C26k7seeLl1CqzwcwmYL0N+C7fEQUenUXgAqh/3XDLRw7mup0Uj17FUOp0gFO29lrIJzt7UmzEBoSjACTCUXDC6FFs2hs2rZHrtTndXzw+VwkJiUjNDgYwcGBUpZFf+s2/4GWzaKxddc+jBz6RLbQi9X0lb/sQGZmJno+8qBdobcey/WEW2hSv4bs157Qz5i/GoXCQhASHISBj7eHwWCQFSEpKWn44++jeKZPJ/y07nd50yGEPi4+ASdPX5A3D+Im5b9/H5XuSqG3omt5Q6zlpVGsoXfr54onaUiAQq8hXC91TaH3EkgvdJORkYxjf85E7IU/snuzCP26sW2R2qA+6jZ/GWGFS3shGrvwFgEKvbdIOu9HhdBv2ZyIwwft189bj9BToZ8+ZwV6d20jpd5oNOLCpatYvXFHLqFfvOJXzFq0FvXrVMF7rz6LGfNXoVPbe2UZjZDpouGFsWzNFoQXKYTa1SvKz0cN740rcdcxeeYPMJmMKBlVDO+MekrePIjDdoXeeizHTp7DuYuxUsjtCf2r7/8HZjPkivvg/l0QGGCSQi+O8MKF8OeB47idlIIBPR+WQi+Oq9cS5Ar+6OG95fgo9KIQyuawvDRK7Es/bcLLaNuikfOfCD9swZIbP0waAAq9/vNGoddPjuIv/Yn9uz7LMSBroY+vHIEK1R5Btej++hk0R8JdbhReAyqE/tCBFGz9LdHprIb+KwoB7i3Q5ypzsUivKEkRK/RC1sUhaurF6r31YVtyI0p1RN39fY3rICgoEHv3H8N9TeqiVFSELGnp2aU1Pp6+EE0b1EKH1k3tCr11+Y+rJTdLV/0m+3q6TydZciOOB5rXl7XxZy/E4v3XnskW+sb1aiDhVqK8ebFUl3CF3sElZnl7rPh43tSxfvnyqLx+eij0Tn+36LIBhV6XackxKAq9fnJ06L/TEHt+V55CHxRSFC07f6mfQXMkFHqF14AKoU9NNWPh/GtITsq1hpo90+gGoWj1YE7Rzg8GeyU3YrV736ETmDZrOWpUrYC4+Oto1rA2RL269WFbctO0QU2cOnsJwwZ2lc3+OXUei1dsRtcOLXD0nzOyhl5I/7hPZ+Pp3p0QXbuKxzX0ovY/MDAw+0YhMPDOyr9YVBY3BCLW9Ikjs4XesthsKRcfMeRxTJ/9kyy7GTtigPxm4W46nG5bKWCIrz1eHDtFvhFWrNhbHrbwd1AUev/MIIVe/3mj0OsnR7s3vIakxEt5Cr34sGXnaQgKKaafgd/lI2HJjboLQIXQi9mcP5eOVT8l2J1YiZIBEOU2RqN28xa7yYSGBMmSGl8dYgv01NS0HOFF7XxB8UpfcRVxcwm9ZVVe1DqJI7xwGL75dIx88GDsxJlYsX67rKOyvbvz5STcjU2hd5ecb8+j0PuWvyvRKfSuUFLTZt+OSbh2ZX+eQm80BeHBrt/AYNDQJtRMt8BEodCrS6UqoRczEvvR//1nMs6eScPNhEyUKRuI2nVDUKny3fFgOoVeu+s6h9BbaucH9eqQXTcvvub4ev5KfDh2WPbON+Mnz8Mrw3r5fQkOhV67C0vLnin0WtL1Tt8Ueu9w9EYvMYeXQOxwY33Y1tCHF6+OJq3f9UY49uElAhR6L4F0oRuVQu/CcNiEBNwikGvbyoK877wtIQq9W9eMz0+i0Ps8BU4HQKF3ikhZg/S0ROz57W2k3I7Ljmkt9NeqRKFx63cRHlFF2ZgYyDkBCr1zRt5qQaH3Fkn240sC+V6h9+VgvR2bQu9tomr6o9Cr4exJFAq9J/S8f27CtRP4c8t7AO48kGct9EUfHY6KNbp5Pyh79IgAhd4jfPk6mUKfL1xsrFMC+aqh1+kc3B4Whd5tdD49kULvU/wuBafQu4RJaSPxQqm4i3twPe4gonb+iei3Z+HGpk1IapD7VexKB8ZgdglQ6NVdGBR6dawZSTsCLu1yo1143/bsb0KfYkpBhiEDZoMZAVkBCMkMgQG+e1rdV9mj0PuKvOtxKfSus/JFy4gtGxD6eHfc3roDCXUa+mIIjOmEAIVe3SVCoVfHmpG0I0Cht8NW/HA7kn3tUuG450xDFq4FXYMQeusjwByAqJQoBJrdfAuFLybjhZgUei9A1LgLCr3GgD3snkLvIUAFp1PoFUD+XwgKvTrWjKQdgWyhF9tVijfDPvpw/r5+Xb1xp9yfvkSk/+1f7A8r9KLiNTbkCtKM6XavApPZhFIpJSH+/245KPT6zzSFXt85otDrOz9idBR6dTmi0KtjzUjaEcgh9INemoCypSJdenmUZYvLi1fi/fYtsv4g9DcDbyEh0P6LKCyXRVhmKCJTI7W7SnTWM4VeZwmxMxwKvb5zRKHXd34o9GrzQ6FXy5vRtCGQq+TG8vKoCmVL2hV16xdPde/YEhPeHKrNyBT06g9CHxd8NVepjS0asTpfNrmMAmL6CEGh10ce8hoFhV7fOaLQ6zs/FHq1+aHQq+XNaNoQcFhDv3nHX3hx7BS7UadNeDn7xVPaDEtNr/4g9BdDLyHTkOkUiBD6u6XshkLv9HLweQMKvc9TkOcAKPT6zg+FXm1+lAr9hdswrzkLxCYDSRlARDBQqxgM7csDIXdP6azaDN8d0fhQrJ086+mhWK7Q504QhV7/v5wo9PrOEYVe3/mh0KvNjxKhT8mEed05YOcVy+sgck6ySCAM3SoCDe6e8lm1WS740Sj0Ohf6W4G3cIM19DmyRKHX/y8mCr2+c0Sh13d+KPRq86NE6C8nwfzZgbwnVi0chmG13Z58UnIqxk+ei0ux8XJL68ji4Rg7YiASbt7Gh9O+x+XYawgJDkTfHu3Qo1OrHHE++WoxDhyJgcloRFhoCEY/1xuV7ymDDVv3YOX67bK8ukjhMByPOY8JX8yH2WxGZEQ4xrzQT36+a+9h3LiZKM8NCgzAkP5dZNve3doiLDRY/vn4yXN4tH3ujVdsx92+dVM81vkBhAQHYdKXC1GuTAkM6PkwUlLTMO7T2WhUtxrubVQb46fMk3N46+VBcqyWwzbW6g07UaNqBRyPOYfvFq1F4+jqGDW8D76cvVzOOSMjE493eVDG/OfUhRzjtJxbo0p53EpMgihL79axJdo/2FSGczS+Pt0fgjXT6lXKY+TQXggNCXI7v85OpNDrXOjNMCMuJA6pxjS7ueQuN84ucX7uCwIUel9Qdz0mhd51Vr5qyV1u1JFXIvSbL8K89pzTSRneawqEuld6IzYrmblgDYYO6CLFWhxCQoX4jhj8OMqVjkJmZhZizl5E9crlc4zli29/RL8e7XLsWJialo4pM39ARmYWOrVthsbRNXDo2GkcPBoDIazbdh/A/sMn8MIzj8m+hOh379gKdWtWgu1YrM+zhWDdNjQkGOs2/4Ej/5zBK8N6YdL0RYi9eh1vvjQA5y/F4au5K9CsYW0M7vcIRGm4ONq2aJSjS9tYi1f8inq1qshxWc/T8mdxozJ+8jwMH9QVN28lZc9PdGp97s49h/D7fw/IG4DXnu+LAJMJH3250OH47DF1egF40IBCr3OhF8PLMmQhIfAmEgMSc4w2JCMExdMj7praecvkuULvwU+8olMp9IpAuxmGQu8mOIWnUejVwVYh9ObFJ4G9V51OyvBiXeCewk7b2WsgxPjtSbMgpFjIZtHwQmjRLBqbtu2RK/V5HR98PheJSckIDQ5GcHAgnuzVUUq5kOuWzaKxddc+jBz6RLbQi5X3lb/sQGZmJno+8qBdobcey/WEW2hSv4bsNy+hFzcit5NS8O8p8zBqeG98v3wjypcpgZTUdFy/cRMloiLkTYozoZ888wd5AyOOC5evyrE7EnqTyYTPvl6CV5/rI9vaO7dOjYqY9t1yPPxAEyxfu03e0FStWFbeIDgan4VppQql0bf7Q4goWsStvLp6EoXeD4TeMsQsQyZSjKnIMpjlW2ID7qK9563TRKF39cfbd+0o9L5j70pkCr0rlHzbhkKvjr8Sof/xFLA71umkPBX66XNWoHfXNlLqjUYjLly6itUbd+QSerHyPGvRWtSvUwXvvfosZsxfhU5t75VlNAaDAUXDC2PZmi0IL1IItatXlJ8Lwb4Sd10Kr8lkRMmoYnhn1FPy5kEctiv01mM5dvIczl2Mld8COBN6cSMhVszFCv3CnzahR6cHMOWbH9CwbjWUL1sCMWcu5RD65o3q4N1PZmH/4Rg827ezXI3f/ddhdG3fQoZatWEHRBtboRfCnXDrNk6evoD3X3sW0bWryBsWe+eWKhEBMR/R/+4/j0B8e9H/sXZS6B2NT7QR5Ud79x+XNwGCVXBQoNNrwN0GFHo/Enp3k1zQzqPQ6z+jFHp954hCr+/8iNFR6NXlSIXQY2cszMtP5T0pA2AY3wwINLo1eXslN1evJWDi1AVS6IWsi0PU1IvVe+vDtjxErIKLuvv7GtdBUFAg9u4/hvua1EWpqAhZktKzS2t8PH0hmjaohQ6t79STe6PkRqzQ7z98EsvX/Y7/e3kgps/+Sd4EJKekomiRwvjz4PFcQu9pyY14tmD2knUYN/ppnL0Qa7fk5sz5K7IGX9zcCM7//fso3ho5CN98v8bh+Cx8L8ddw7RZy2UOxPMEWh0Uegq9VteWZv1S6DVD67WOKfReQ6lJRxR6TbB6tVMKvVdx5tmZEqFPzoT5k33ALftvfZcDbFX6zk43bh72Sm4G9++CfYdOSKEUD4bGxV/PrkG3DmNbctO0QU2cOnsJwwZ2lc3+OXUei1dsRtcOLXD0nzOy5ERIv3hI9eneneTqtidCL8pzggIDcersRVSpWBavv9Bf3nTY3miIunmxQi++Tfjoy+/l2ERbS3mN+LujGvprN27KbxfatWoMweXreSuznxv4Zcse7Nl3FF3a3Y+jJ+7MTxzim4xa1SvK0iOxIi/e0SQeCP7s66Vo90BjWYpk/eyBZXzioeN3P/4OEcWKIObMRTzdp1P2g7RuptfpaQ6FXjzJK+5axJ7z4hB70u/+6wjCC4fhm0/HyK8u/P3wh33o/Z2xFuOn0GtB1bt9Uui9y9PbvVHovU3U+/1R6L3P1FGPSoReWnECzN8ctb9tZflCMLxQFzAZNJu42E1G7LIiSmp8ddxMTELq/9utxvooFBaS/RCvr8alVVxx0xMSEozAAPcedM7PuOwKveVtsK+/2E8+PSzuOOYt/UXKvZB6y58tT1HnJ6Ce2lLo9ZQN18dCoXedla9aUuh9Rd61uBR61zj5shWFXh19ZUIvpiT2o99yCTh2A4hPASqHw3BvCaBOhLoJ+zDS3Sb0KlE7FPrn35ws64nESrxYrReH2IdUfJUhvmKZPnFkju2NVA7aW7Eo9N4iqbYfCr1a3u5Eo9C7Q03dORR6dazdjUShd5dc/s9TKvT5Hx7PIAGXCNgVelGHJUpsBvXqgHo1K2PQSxNgvVr/0bSFmDd1LIXeJcRs5G0CFHpvE/V+fxR67zP1Zo8Uem/S1KYvCr02XO31SqFXx5qRtCPgsIZerMQPGT0J4uuR7h1bytV5SymOeMuW+Lu/H1yh988MUuj1nzcKvb5zRKHXd37E6Cj06nJEoVfHmpG0I8BdbuywFT/cjmRfu1SwZ1cJUOhdJeW7dhR637F3JTKF3hVKvm1DoVfHn0KvjjUjaUeAQk+h1+7q0qhnCr1GYL3YLYXeizA16IpCrwFUL3dJofcy0Dy6o9CrY81I2hGg0FPotbu6NOqZQq8RWC92S6H3IkwNuqLQawDVy11S6L0MlEKvDigj+YSASzX0tiMTG+vzoVif5ItBAVDo9X8ZUOj1nSMKvb7zI0ZHoVeXI67Qq2PNSNoRyHOXm5b3RsvX/n49fyU+HDtMbvwvtrBs37qp3J/e3w8+FOufGaTQ6z9vFHp954hCr+/8UOjV5odCr5Y3o2lDwOk+9CKs9b7z1i+Z4oultEkKe82bAIVe/1cIhV7fOaLQ6zs/FHq1+aHQq+XNaNoQcCr0JaOK4fXxX+Ojt4bLfef5YiltEsFeXSdAoXedla9aUuh9Rd61uBR61zj5shVLbtTRVyn05+IysHLnbVy+nonbKVkoXsSEuhWD0LlZGEKCDOomzUgFjoDTkpvB/R6RZTZVK5WD+PO3C3/G9j8OYNqEl2UJjj8fLLnxz+xR6PWfNwq9vnNEodd3frhCrzY/KoQ+Jc2M1btvY9uBFJjtTC88zIjHWxVCo2rBaifPaAWGgEu73FheKHXuYizCC4fhm0/HoG7NSn4PgULvnymk0Os/bxR6feeIQq/v/FDo1eZHhdBfjM/Eh4uv5zmxmuUD8UK3om5PPik5FeMnz8Wl2HgYYEBk8XCMHTEQCTdv48Np3+Ny7DWEBAeib4926NGpVY44n3y1GAeOxMBkNMrF2tHP9Uble8pgw9Y9WLl+u3yZaJHCYTgecx4TvpgPs9mMyIhwjHmhn/x8197DuHEzUZ4bFBiAIf27yLa9u7VFWGiw/PPxk+fwaPv7c83PetwZGVl47sluaNmsHnbuOYQvZi1DaHCQPKdnlwdRo0oFfPHNj3jrlUEoXaI4dv15GDPmrcLH7/wL4UUK4fvlG7Fi3e+yffdOrdD/sYdx4VIc3vl4Fm4kJKJhvWp47V99ERISjMkzl2Lbrv2oX6cqXn+hn5xfSmoa5v3wC3o+8qCcnzis2VSvUh4jh/ZCaMidMenpcEno9TRgb46FQu9Nmur6otCrY+1uJAq9u+TUnEehV8PZkygsufGEXv7OVSH0G/5Mxqpdt50O7KPBkQgNdq/0Jik5BTMXrMHQAV2yKyhuJSZh/JR5GDH4cZQrHYXMzCzEnL2I6pXL5xjLF9/+iH492snSasuRmpaOKTN/QEZmFjq1bYbG0TVk2fXBozHo0/0hbNt9APsPn8ALzzwmTxGi371jK7ngazsW6/NsIVi3FXL/6VdL8MZL/fHngeOyqfUmLKKfad8tx0OtGuGxzg/ItsdjzsmNW/779zFcvHJVVpOIQ1SUlC0VhY5tmsFoNMBgMGDmgtWoVe0eBJhM+OfUeTzZqyNWrN+OwIAAtHugMb5fthG///eA7M/Cwh4bp4n0QQMKvR3ofFOsD67EfISk0OcDlo+aUuh9BN7FsBR6F0H5sBmFXh18FUI/f9Mt/HEs1emkRj1eDJVKBThtZ6+BEOO3J81CaEiwFNai4YXQolk0Nm3bI1fq8zo++HwuEpOSERocjODgQCm6or91m/9Ay2bR2LprH0YOfSJb6MXK+8pfdiAzM1OuZtsTeuuxXE+4hSb1a8h+8xL6y3HXpVS/9nxf7NhzEAuXb5JSLg4h3MWLheO/+47i0pV4NG1QE1firuPw8dN4ZVgvKeuWGwrRXsj/mo075bcI4kjPyMTkGUvRo3Mr7NhzCPeUKylvFsSNw5F/zmJAz4dlO1uBt7CpVKE0+nZ/CBFFi7iVH61PotBT6LW+xrzeP4Xe60i93iGF3utIvdohhd6rODXpjEKvCVa7naoQ+sVbErH9UIrTSXkq9NPnrEDvrm2k1BuNRly4dBWrN+7IJfSLV/yKWYvWon6dKnjv1WcxY/4qdGp7rywzESvZRcMLY9maLbKMpXb1ivLzUcN7S4GePPMHmExGiE1T3hn1lLx5sCf01mM5dvIcRNm2+BbAntAL+Rer6GfPx2LyBy+hTMniELsq3ryVhBZN68pTCoWF4NTZy/IbgnJlSmDWwp/xfy8PlCvxjoR+xfrf5dxFidCcpethNBgw6IkOcu5VKpaRQn/g6CnZp2VstkIvvqkQZUR79x/H8rXb5JyDgwKd5lJ1gxxCb6mVf6ZPJ3y3eJ2Eb+/gi6VUp4nxrAlQ6PV/PVDo9Z0jCr2+8yNGR6FXlyMVQv/7oRQs2ZKY56REoc3Hw6IQ5N4Cfa4yFxHs6rUETJy6QEqtpSZc1NSL1Xvrw1ZiRamOqLsX7yIKCgrE3v3HcF+TuigVFSHlt2eX1vh4+kI0bVALHVo3tSv01uU/rpbc/P7HQezZdxSvPd8Pv/+xX/ZrW3Ij4ov6+DPnr6BC2RIYP3meFPqN2/aicFgounZoIc9b9csOpKSloUenBzB78VqUjIpAtw4t5A3Ltt37cf5SnJR48ZzA9Ru3ZL2/OByV2FyOu4Zps5ZLluK5AL0dXKG3kxGW3OjtMs05Hgq9vvMjRkeh13eOKPT6zg+FXm1+VAh9UqoZ/154DbeS7O1xc2e+reuHyp1u3D3sldwM7t8F+w6dkCJao2oFxMVfR7OGtbPrzC2xbEtuRDnLqbOXMGxgV9lE1JsvXrFZyvLRf87IGnoh/eI9RU/37oTo2lW8UkMvvlmYvXidjFmxQiks+unXXCU3lhp+0UbM2SL0IcFBmDR9EZJT7pQ2ib7GPN8Xe/Yfw5ff/YR6NSvLfxelOw3rVsO/v5gvV+yv3biFca8+LUtpRLnP6o07Ua9WZTl30ee7H3+HiGJFEHPmIp7u0wntH7xzA6O3I0+hF19jfPb1kuwxF6QdbsSk+FCs3i5H18ZDoXeNky9bUeh9Sd95bAq9c0a+bsEVenUZUCH0YjbHzqdj+soEu9tW3lMyAK/0LAaTUbt5iwdOxe4sYoXaV8fNxCSkpqblCC9Kaby5DbqQfJPJ5LQsRpThCCZitd0ZE3HzInbGCQy4U16kx8Oh0AuZX7rqN8ybOjb7SV/xlcmQ0ZMwYezQHF+B6HFiroyJQu8KJf21odDrLye2I6LQ6ztHFHp950eMjkKvLkeqhF7MSOxHv+nvZBw5m4a4hExULROIFnVCUK+S/rZB1CIDKoRei3H7Q58O3xQ76KUJeP3FfrnEXTykMG/pL3yxlD9kt4COkUKv/8RS6PWdIwq9vvNDoVebH5VCr3ZmjHY3EXAo9M+/ORnjRj+d6wVSYpVe1ExNnzgyx36lvoYmvmJ5cewUOQzrt9halw01b1Q7x2dcofd11tyLT6F3j5vKsyj0KmnnPxaFPv/MVJ/BFXp1xCn06lgzknYE7Aq9RY4H9eqQa4Vej0JvGe/uv47AWtrFtwkfTVuYXTY0duJMSVK88UwcFHrtLiwte6bQa0nXO31T6L3DUateKPRakfVevxR677F01hOF3hkhfu4PBBzW0DsqrREr3idPX8iWYj1MUoh61Url5FC2/3EgexXe8u+Wt4bZCj6FXg/Zy/8YKPT5Z6b6DAq9auL5i0ehzx8vX7Sm0KujTqFXx5qRtCOQLfSWPegd7T1vPQQ97UNvveoubjYsQi/GK0pwWt4bnb09k/h2YdS4L/HZuBdkKdGV6/Zf8lAqIsThZ9qlgj27SiA40IjQIBNu3E539RS2U0ygWKFAJKdlIjU9S3Hk/IUzwwwDfLfjQ/5G673WRX/7BSE9u+P2tp1IrNvQex2zJ68RCA02IdBkwM2kDK/1yY7sExD/zedBAv5OwK/3obf9tsCe0FuXDdkKfUamfdkIMBnh6DN/T3hBGL/YXkrsupWV5Xg/34IwT3+eg3jjn9kM+XY+PR9pGVkICtBwnzidTt649mcYu3ZF1q7dyGqqzz2VdYpO2bDk7zkAWTr/GVIGRMNA4r/5PEjA3wn4tdCL1fkV67fnyoGoo5/09nCM+eDrPFfoWXLjn5cvS270nzeW3Og7Ryy50Xd+xOhYcqMuRyy5UceakbQj4NdCb4vFeoVevKSANfTaXTi+7JlC70v6rsWm0LvGyVetKPS+Iu96XAq966w8bUmh95Qgz9cDgQIt9NzlRg+XmPfHQKH3PlNv90ih9zZR7/ZHofcuTy16o9BrQdV+nxR6dawZSTsCBVroBTbuQ6/dxeOrnin0viLvelwKveusfNGSQu8L6vmLSaHPHy9PWqsU+lvm0ziRsQi3zReRjkSEGCIRZWyEKqbHYEKoJ9PguXc5gQIl9PnNJWvo80tMH+0p9PrIQ16joNDrO0cUen3nR4yOQq8uRyqEPhPJOJGxFOezNgLIvVlAEIqhRsBAlDI2VzdxRipQBCj0dtIpfrgdyX6Byr6fToZCr//EUej1nSMKvb7zQ6FXmx8VQp9oPofd6WPznFhxQ100CnzD7cknJadiycrN6B/bn0sAACAASURBVN2tLcJCg3E85jyOnzyHR9vfn2efopKhSsUyaHN/Q4gXdM6cvxpiK/OIYkXQ69E2eKhVY4yfPBeXYuPlNr+RxcMxcmgvTJ21DJ3bNkfr+xsg/vpNjBn/FYb064IypSLxzsezcCMhEQ3rVcNr/+qLIoXD7I7h9LnL+HDa97gcew0hwYHo26MdSkVFyLb3N60r/1/M6/3P59iNZWnjNrQCdKJDoc9rX3o97UPvSS64Qu8JPd+dS6H3HXtXI1PoXSXlm3YUet9wz09UrtDnh5ZnbVUI/enMVTiZucTpQFsHfY0A2JdfZycnJadg5oI1GDqgC8TGIGKr7oNHY6Tg/33oBGKv3kCJyKJoUKcaTCajFP5TZy9h555DaN2iAUpGRsidA0cN74WQ4CBkZmYh5uxFlCsdlaPfO5KdgolTv5dbSI8dMRBbdu7D3KXrMWRAFzzYvAHE1sVi69WZC1ajVrV78EDz+rmGfysxCeOnzMOIwY/LGJZ45y/FybZtWzT6n9A7jiXanL1wBUf+OYvyZUrIWNdu3MT+wzGIKFZYzvXqtQRcjrsm/71YeGGUK10CB47EoEK5kqhRpbwzrH7zuUOht35hk9/MJp8DpdDnE5hOmlPodZKIPIZBodd3jij0+s6PGB2FXl2OVAj94YyvcSnrd6eTahr4LooaqjltZ6+BkOy3J81CaEgwAkwmXE+4hSb1a6B/z4cRc+YSIooWxg9rtqBiudIoVrQQft60G8MGdsWcJevQqnm0bCNW6i0ibYlh22/R8ELo91g7LF6xGeGFw1CqRHGcOH0eJqMRdWpWyj4/PSMTk2csRY/OrVC9cm5xFjccK9b/Lm8IrA+xoYmt0IsbFXuxgoMCsWnbn/IFooeOn0bViuUwe8laPNW7E/b8fRQBAQEoHlEES1b+hleG9ZJzFTcrgsm336/BqOG9EVW8qFu89XaSXaEXq/PPvzkZ40Y/Ld+oWlAPCr1/ZpZCr/+8Uej1nSMKvb7zQ6FXmx8VQn80YxYuZG12OjFPhX76nBXo3bWNlPpjJ8/h3MVYPPFoGyxd9Rs2bdsrS2O6dmghS3Lq1aqC6FqV5eYhQuSthX7SlwuxdvNutGwWjbEjBsC6X6PRiKDAAMxatBaPdW6FcZ/MhniJp/X54qWCc5auh9FgwKAnOsjVetsjv0JvL9aJUxdwX5O6ch7isHwr0af7Q7Js6Ku5K9G8SW2cuxAnpd8yV3HTMuGL+ejesVWB8VwKvZ0fL9bQO/2d49MGFHqf4ncpOIXeJUw+a0Sh9xl6lwNzhd5lVB43VCH0F7I24WjGbCdjNaBt0EwYEezWnByV3FSvUh679h7Gv57qjt92/i3FOzDAJFfNRQ26RXIDTAHYuecgRj/XR5bkWIR75NAn7JbcWMp7zl2MQ8XypbBg2UZ5Y9Dq3vqYvXgtSkZFoFuHFnZlXkxQlMJMnLpArtBHRoTLOSfcvI0/Dx6Xf7YuuXEU69yFWFS+pyweaB4tzzl64ix27DmIZ/s+Im9mFv30K5o0qIFTZy/fnUIvoNi+lMmtq0vnJ3GFXucJcjA8Cr3+80ah13eOKPT6zo8YHYVeXY5UCH0GbmNn2utIQ4LDiVUwdUQNU87yk/xQcCT0rVs0xNsfzUKpEhE4c/4K2rRoiJbN6uH9z+agbOko+eDsK8N7SRFfsGwDNm//C5UrlMHluHg0a1gb/Xo8lKOUx7rkxlKvL8ZpuTEQf/7yu59Qr+adVfN2DzS2W0MvPtu2ez+mzVqOGlUrIC7+uownbgoWLt+EsqWi5PndOrbAtt0Hsp8NsI5Vs0oF+cBsicgIpKamYdRzfTBv6XpcuXodaenpsqRIrNSLm5i7coVewBJ3Zl/PX4kPxw6TD1cUxINC759ZpdDrP28Uen3niEKv7/xQ6NXmR4XQixldMx/EX+mT7G5bGW6oDFFuY4BJk8mLevb09AxZamM57P2b+EyUy4idZURbe6Uy7gzwZmKSFG7ro1BYSLZfinihIUFeiyceuA0JCZbfRNwth8OSm0EvTZBfV9g7uMvN3XJ56HOeFHp95sV6VBR6feeIQq/v/FDo1eZHldCLWYn96M9k/oz4rH1IMl9BMWNNlDO2lS+XKsiHM6EvyHNXNTfuQ2+HNGvoVV1+7sWh0LvHTeVZFHqVtPMfi0Kff2aqz2DJjTriKoVe3awY6W4jQKGn0PvdNU+h13/KKPT6zhGFXt/54Qq92vxQ6NXyZjRtCOQp9GIv0BfHTskRedqEl3PtUarN0LTvlTX02jPWIgKFXguq3u2TQu9dnt7ujULvbaLe748r9N5n6qhHCr061oykHQGHQi9k/qNpCzFv6liUiCwmRyAelB0yehKGDHhUPi3s7weF3j8zSKHXf94o9PrOEYVe3/nhCr3a/FDo1fJmNG0I2BV6sfWRWJkXLwqwfWOYEP15S3+BWKn3991vKPTaXFRa90qh15qw5/07E/oscxri0pchKfMwMrJuwGgoglBTRZQKGgCToZDnA2APeRKg0Ov/AuEKvbocUejVsWYk7Qjk+8VSYpV+3KezMX3iyOyVe+2Gp23PFHpt+WrVO4VeK7Le6zcvoU/LisXF1BlIzbqQK2CAIQJlQ4Yh1Fhw31DtPcru90Shd5+dqjMp9KpIAxR6dawZSTsCXKG3w5a73Gh3wXmjZwq9Nyhq20deQn8h5T9IzNzvcADBxrKoGPomDAjQdpB3ce8Uev0nn0KvLkcUenWsGUk7Ag5r6MUbv5au+o019NqxZ89uEqDQuwlO4WmOhD7NfBmnkt5zOpLyIS+ikKmu03Zs4B4BCr173FSeRaFXR5tCr441I2lHgLvccIVeu6tLo54p9BqB9WK3joT+VsYeXEz91mmkqKDuiAzs5LQdG7hHgELvHjeVZ1Ho1dGm0KtjzUjaEeA+9BR67a4ujXqm0GsE1ovdUui9CFODrij0GkD1cpcUei8DzaM7pUIfvx/mveOBhH+A1OtAofJA+YdhaDgaCCyibtKMVOAIUOgp9H53UVPo9Z8yrUtuAk8cQsDJQzBev4rMqNJIr9UImeUr6x+MTkZIoddJIvIYBoVeXY6UCH36LZj/nAgcnQ3AnHtyoSVhuPcDoFI3dRNnpAJFIIfQx8XfwKCXJuCZPp3w3eJ1OHcx1u5kK5QtmaO23l+JcJcb/8wchV7/ecvroVixw82tjL8cTiLIWAaVQv8PBphytTEmXEORbz5E6PolOT8zmXD78aFIHDgC5pAw/QPy8Qgp9D5OgAvhKfQuQPJSEyVCf/0IzCsfynvEZR6AoYPN77Z8zDEpORVLVm5G725tERYajOMx53H85Dk82v7+PHsRz0xWqVgGbe5viN1/HcHM+ashfDCiWBH0erQNHmrVGOMnz8Wl2HgYYEBk8XCMHNoLU2ctQ+e2zdH6/gaIv34TY8Z/hSH9uqBMqUi88/Es3EhIRMN61fDav/qiSOHcv5fNZjOWr92GxSs2Izk1DWVLReLtV57E2fNXsHT1Fox//Vm5Pfqyn7di595DeO/VZ7F4xa/o1rElIiPCZcyV67ejR+dWeO/TObh45SreeLE/GkfXyDHfT75ajANHYpCRkYnHuzyIxzo/ID9ft/kPzP3hFyQnpyIqsijefGkAqlYsiw1b98h+J7w5NMe4Dx47Jfm+8eIAyVdvB1fo7WSEu9zo7TLNOR4Kvb7zI0aXl9Bnmm/jdPIHyDAn5JqIwRCEiiGvQ+x0Y+8o9u5QhOzc4BBA4qCREP/jkTcBCr3+rxAKvbocKRH6A1Nh/nOC00kZ+h0Fgoo6bWevgXiH0MwFazB0QBcpwmKb8YNHY6Tg/33oBGKv3kCJyKJoUKcaTCajFP5TZy9h555DaN2iAUpGRmDF+u0YNbwXQoKDkJmZhZizF1GudFSOfkVsEWvi1O9hMABjRwzElp37MHfpegwZ0AUPNm8Ao9EAg8GAmQtWo1a1e/BA8/q5hnzyzEV5zlsjn0RggAnJKWkwm7PkTcXC5ZvwZK+OaFCnKiZOXYDklFT8+40h+Ob7NejXo53cNl3cdCz8aRNGDH48e659uue+afri2x/lOeKmYvzkeRg+qKu8GVizcSfGvNBfxk64eRvJqamIKFoEU2b+gIzMLHRq2yz75uDi5av4cvZPksvo53rr8j1M+d6Hni+WcuvnjCd5kQCF3oswNerK2YulxFfON9J/R2Lm30jOPIVgUzkUNkWjWGBbGBFod1RBh/ei+MjH8x6xwYDYhX8gq3gJjWZWMLql0Os/jxR6dTlSIfTm318GTjpffTc8shoo0cStyQvJfnvSLISGBCPAZML1hFtoUr8G+vd8GDFnLiGiaGH8sGYLKpYrjWJFC+HnTbsxbGBXzFmyDq2aR8s2YqXe9oWitv0WDS+Efo+1kyvr4YXDUKpEcZw4fR4moxF1albKPj89IxOTZyyVK+jVK5fPNSexgj/y3Wlo06IhWjWLRsUKpaVcC8+Mi0/AydMXULdmJXlj8d+/j+KtkYPw6VdLkJiUjNDgYCnghcNC5aq+5eYlL6E3mUz47OslePW5Pvjx56125/rPqfNy5b5ls2hs3bUPI4c+gVuJSZi/bINc2RdzttwwuZUkDU/Kt9DzxVIaZoNdu0SAQu8SJp82ci70+R9e6OoFKPrF/zk98drHC5HWIO+vmJ12UsAbUOj1n2AKvbocKRH6nWOA4/OcTspToZ8+ZwV6d20jpf7YyXOydPqJR9vIbcg3bdsrV6a7dmghS0bq1aqC6FqVYSm5sRb6SV8uxNrNu6XYjh0xANb9Go1GBAUGYNaitXiscyuM+2Q2BvXqkOOGQJTTzFm6HkaDAYOe6CBX6+0dqWnp2HfoJHbsOShX5se/PhhnL1yRTa9eS8DWXfsxenhvzJi/Sgr9jPmr0antvdklN0K+hXRbC70oyxFjq1+niizT+fSrxUi4dVveILz/2rOIrl0le862Ny/i3PAihVC7ekUZc9Tw3li9cSca1q2G4sWKYMGyjXimb2eULlHcaS5VN8i30IvEb//jAKZNeFmXXznkByBr6PNDSz9tKfT6yYWjkWgh9EUnj0Xoz987nXzCiH8j+dEBTtvdzQ0o9PrPPoVeXY5UCD2OzYF51xtOJmWAYcBJICDUrck7KrmpXqU8du09jH891R2/7fxbirdYCRer5vc3rZsttwGmAOzccxCjn+sjS3KEJK9Y/7sUZutSHjE461jnLsahYvlSUnbFCn+re+tj9uK1KBkVgW4dWjiUeVFiYzIaEBR051tZS3lORmam/HvjejWQcCtR3px8PmOpFHpPSm4ux17D7CXrMG700/j70Ens/uuwvFkQNxvi24Rr1xPwxbfLcF/jOnJMe/cfw31N6iIqIhyXYq8hLS0dG7ftxcDH26NZw1pu5UjLk3IIvUjekNGTcDMxyWFM8fXKN5+OkV+D+PtBoffPDFLo9Z83LYTe1RX6+CnLkF67sf4h+XCEFHofwncxNIXeRVBeaKZE6NNuwPzTg0BynOMR1x5yZ6cbNw9HQt+6RUO8/dEslCoRgTPnr8gSl5bN6uH9z+agbOko+eDsK8N7SRFfsGwDNm//C5UrlMHluHg0a1gb/Xo8lKOUx7rkxrr8xLLSL4b/5Xc/oV7NOzuPtXugsd0a+tPnLssHae8pV0q2E+I+5vm+2LP/mPy7ZfVc1MrnJfRP9+mMT/6zCKfOXpbfJojVdevDUkMv6u5/2bIHe/YdxchhvTFj3krsP3wS5cqUkN9kPPrw/bh246YsQxKHKL8RJTZvvNRfljDZ8nUzTZqdlu8Ves1G4oOOKfQ+gO6FkBR6L0DUuAsthD7gXAyihncAMjIcjt5cqAhiF+yCOayQxjP07+4p9PrPH4VeXY6UCL2YzqWtMP/S1/62lZENYXhkJWC0/wyRpzTECnR6ekaO3Vns/ZuII8plxI45oizHUalMfscjFopTU9NynFYoLERWeoiyG3EE/2+lPr99e9Je1OenpWcgNCTIk250cS53ubGTBu5yo4tr0+EgKPT6zo8YnRZCL/ottGg6isya5BDAjXEzkNKig/4B+XiEFHofJ8CF8BR6FyB5qYkyoRfjFfvRH/wPcGEzcOsUUKo5DNUHABUK9u+tvITeS2m867uh0FPo/e6HgEKv/5RpJfQwmxH6yw8oMnMCjDevZ4PILHsPEkZORFrDlvqHo4MRUuh1kAQnQ6DQq8uRUqFXNy1GussIOBT6vOrp+WKpu+wq0dl0KfQ6S4id4Wgm9JZYqSkIOrIXATFHkVGzIdJqNgACAvQPRicjpNDrJBF5DINCry5HFHp1rBlJOwJ2hV4U/r84dgpa3hstn/b9ev5KfDh2mKx1GjtxJtq3bpprn1Lthqhdz6yh146tlj1T6LWk652+NRd67wzzru2FQq//1FPo1eWIQq+ONSNpR8DpQ7Ei9LhPZ2P6xJHyzVx8sZR2yWDPrhGg0LvGyZetKPS+pO88NoXeOSNft6DQq8sAhV4da0bSjoBToS8ZVQyvj/8aH701XAo9XyylXTLYs2sEKPSucfJlKwq9L+k7j02hd87I1y0o9OoyQKFXx5qRtCPgtORmcL9HZJlN1UrlIP7MF0tplwz27BoBCr1rnHzZikLvS/rOY1PonTPydQsKvboMUOjVsWYk7Qi4tMuN2NR/0EsT5Mb7fLGUdslgz64RoNC7xsmXrSj0vqTvPDaF3jkjX7eg0KvLAIVeHWtG0o6AS0KvXXjf9syHYn3L393oFHp3yak7j0KvjrU7kSj07lBTew6FXh1vCr061oykHQGnNfR1a1bKEZ0PxWqXDPbsGgEKvWucfNmKQu9L+s5jU+idM/J1Cwq9ugxQ6NWxZiTtCORb6PlQrHbJYM+uEaDQu8bJl60o9L6k7zw2hd45I1+3oNCry4BKoT8Zfxpz/1iMczcu4lZqIkoUjkTTCg3Rt/FjCAsMVTdpRipwBPIt9HwotsBdA343IQq9/lNGodd3jij0+s6PGB2FXl2OVAh9Unoy5u9ZirWHN8EMc67JRYQWxeD7B6JV5ebqJs5IBYpADqHP6+2wllnzodgClX+/nAyFXv9po9DrO0cUen3nh0KvNj8qhP7M9XN4edn/5TmxBmXr4r3Or7s9+aTkVCxZuRm9u7VFWGgwjsecx/GT5/Bo+/vz7FMs1FapWAZt7m+I3X8dwcz5qyE2Q4koVgS9Hm2Dh1o1xvjJc3EpNh4GGBBZPBwjh/bC1FnL0Lltc7S+vwHir9/EmPFfYUi/LihTKhLvfDwLNxIS0bBeNbz2r74oUjgs1xisxxscFIhFK35F1UplYc4yY8aCVZg4dhhKlyiO7xatRbeOLeVclq7egvGvPytfdLrs563YufeQHMvH0xfh4pWreOPF/mgcXcNthv58Yr5X6P15srZj50Ox/plNCr3+80ah13eOKPT6zg+FXm1+VAj9j/tWYd6epU4ntmDQVygUlFt+nZ4IICk5BTMXrMHQAV2k8IpF2oNHY6Tg/33oBGKv3kCJyKJoUKcaTCajFP5TZy9h555DaN2iAUpGRmDF+u0YNbwXQoKDkJmZhZizF1GudFSOfsVYRKyJU7+HwQCMHTEQW3buw9yl6zFkQBc82LwBjEYDDAYDZi5YjVrV7sEDzevbEfo74xVboi9fuw0BASb07f4Qftv5N7bu2o8q95TBoCc64Itvf0S/Hu1w8NgpLFy+CU/26ogGdapi4tQFSE5Jxb/fGIJTZy/Lufbp/pArqApkG+5yYyet4ofbkewXyKvAzyZFodd/wij0+s4RhV7f+aHQq82PCqGfsnUGNv/zu9OJfdT1HdQsWc1pO3sNhGS/PWkWQkOCEWAy4XrCLTSpXwP9ez6MmDOXEFG0MH5YswUVy5VGsaKF8POm3Rg2sCvmLFmHVs2jZRuxUt+2RaMc3dv2WzS8EPo91g6LV2yWW5mXKlEcJ06fh8loRJ2albLPT8/IxOQZS9GjcytUr1zertDPmL8ahcJC5A3EwMfby5sAsflKSkoa/vj7KJ7p0wk/rfs9W+jj4hNw8vQFiA1bxA3Hf/8+irdGDqLQA6DQU+jd+sXhy5Mo9L6k71psCr1rnHzVikLvK/Kux2UNveusPG2pQuj/s/07rD+62elQPRX66XNWoHfXNlLqj508J98f9MSjbbB01W/YtG2vLI3p2qGFLMmpV6sKomtVli8MFSJvLfSTvlyItZt3o2WzaIwdMQDW/RqNRgQFBmDWorV4rHMrjPtkNgb16pDjfLPZjDlL18NoMMhVdiHqtoe4UXj1/f/AbIZccR/cvwsCA0xS6MURXrgQ/jxwHLeTUjCg58NyhV4cV68lyBX80cN7Y8b8VRT6/4Gl0FPonf6C0VsDCr3eMpJ7PBR6feeIQq/v/IjRUejV5UiF0K878iu+2jE7z0mJ+vRFT81AcECwW5N3VHJTvUp57Np7GP96qrssZxHiLsRZrJrf37RuttAHmAKwc89BjH6ujyzJESU7K9b/jpFDn7BbcmMp7zl3MQ4Vy5fCgmUb5Y1Bq3vrY/bitSgZFYFuHVrYlXkxQct4RcmNuOEQx9N9OskxikOU6Yja+LMXYvH+a89kC33jejWQcCtR3rR8PmMphd6e0FveCCu+4vhu8Tp5Z2fvqFC2JOZNHYsSkcXcuuj0chJr6PWSifyNg0KfP16+aE2h9wV112NS6F1n5auWFHp15FUIfWLqbbz44xu4kZzgcGKP1u2AIfcNdHvijoS+dYuGePujWShVIgJnzl9BmxYN0bJZPbz/2RyULR0lHzZ9ZXgvKeILlm3A5u1/oXKFMrgcF49mDWujX4+HcpTyWJfcWOr1xaAtK/3iz19+9xPq1aws59LugcZ51tCLPgIDA/Hx9IVo2qAWAgNN8jxR+mO9Vbplhd5SEiScVQj9iCGPY/rsn2TZjfg2oXb1im4z9OcTuUJvJ3usodf3JU2h13d+xOgo9PrOEYVe3/kRo6PQq8uRCqEXs9l38RDGrZ1kd9vKalGV8WHXdxBgvCOz3j5EPXt6eoYstbEc9v5NfCbKZcQONKKtvVIZd8Z2MzEJqalpOU4VtfPi4V0e3iGQp9CLu63Pvl6SHakgbVkpJsUVeu9cRKp7odCrJp7/eBT6/DNTeQaFXiVt92JR6N3j5s5ZqoRejE3sR//T/p/x5/n9uHTzCuqUrokONdug2T05H0R1Zx56PodCr312HAq9kHlR02RdWmPZp37C2KG5noLWfqjej0Ch9z5TFT1S6FVQ9iwGhd4zflqfTaHXmrDn/VPoPWfoag8qhd7VMbEdCeSXgMN96Ae9NAGvv9gvl7iLp4/nLf0F0ya87PdflVDo83u56KM9hV4fechrFBR6feeIQq/v/IjRUejV5YhCr441I2lHIN8vlrJ+QIEPxbqXmPR0M5JuZcmTw4oYERiYezsn93q+O86i0Os/zxR6feeIQq/v/FDo1eaHQq+WN6NpQ8Cu0IsnpV8cO0XuK2r7ggEKvfuJyMoEDv2VipjjGTk6qVorEHUbBck3rvFwToBC75yRr1tQ6H2dgbzjU+j1nR8Kvdr8UOjV8mY0bQg4rKF3VFojauvFW7omvDlUmxEp7FVlyU1Kchb+2JqK6/F3VuZtj6hSJjRpEYSQUKNCAv4ZikKv/7xR6PWdIwq9vvNDoVebHwq9Wt6Mpg2BPGvoHe1Dbz0Uf96TXqXQnzqegf17UvPMYv2mwahcI0CbTBegXin0+k8mhV7fOaLQ6zs/FHq1+aHQq+XNaNoQ8Ot96C2lQbv/OpJNRzysa10mZL31ZvNGtXM8zKtS6LeuT3a4Om8ZfESkEQ92DNUm0wWoVwq9/pNJodd3jij0+s4PhV5tfij0ankzmjYE/FroLW8Je2vkILnjjigTGjthJr75dAzq1qwk//7RtIXZW2+OnThTUrSUC6kU+p+XJkE8DJvXIR6OfaRXmDaZLkC9Uuj1n0wKvb5zRKHXd34o9GrzQ6FXy5vRtCHg10Jvi0QIvvV2m0Lgq1Yqh8H9HpFNbQVfpdBvWpWMxP/tbOMolYWLGNGuK1fonV3qFHpnhHz/OYXe9znIawQUen3nh0KvNj8UerW8GU0bAjmE3iLEz/TphO8Wr4OjGnq91s2LHXhGjfsSn417AZXvKS136ml5b3S20P9/7J13YFTF9se/27LpBQglAQKh9yKd0LuAgEgvNooFAcWK+OQpD9RnoSpKBwGpEekIAtIF6R0SQkICpDdSdrO7v98Mb2PKbnazuffmJjnzD7A7c86Z75mEz5177tyc37MdfCmB/u8TmXhwP/fpNnlTWjVAjWc6/vNaZnFSXvKtEtDLP4cE9PLOEQG9vPNDQC9tfqQE+qwYPVJPJ8CQkAVjhgEqDzWcAlzg1soLCic66k7azJcub6Vmh95cT28GeEtHb+YF+iyD5RIYtUoBa985mv7kZCO2/ZKALCtM7+SkwPPDveHhSafc2NKYHe+pVChgMBZcwmTLDn0vngIqpQJGkwkmmadIbzBCoyp7P3PKvXugHDgAxtNnYGzVWryFQJYdVkCpVIDhHf2ec1hCuwey//PFbiadCU/OJCL9aqpFV0pXJdw7+kBbm8puxc5FabVfKoDeDO+VK5bLro/PC/gsgXmBPjoxw2JeK3o7w9p3RVkIYSFZ+PukZZ/tuzrDrxqdcGOPvk4aFVyclEh6orenO/WxQ4Fj947hj7t/IDYtFhqlBhU9KmJ40+EILBdox+j8XbzcNEjXGaHTGxwaL9Ugk8kERRl8AYTXkQPQDhmEtGOnkNq4uVRyk59CKMDuRGpUCqSkF3xntxAmqasVBdj/+WI3Q5we8ZsfFehGU1UL74EVHQ4lLT0Tm387jOHPdYOrixa3Qx/gdkgEBvRqX6BNdnhIYEAVdG3fHOyQkWU/7wKr2PDx9sCwAV3RPagl5sxfi4fRcWCXmeXLeWL6xGFYtHI7+nVriy7tmyEuIRnvz1mKQltQKQAAIABJREFUCaP6o0ql8vjXf1ciMSkVzRvXxnuvj4SHe/4LFRZvTru9urTCkH6d4Kx1wldLNsK/ii/GPN8TGZk6zP5mNVo0qo02LRpgzoJ1fD6zpo1DzepVsueWd767fj+FurWq4XZoBFb9shctm9TBO5NHYMnqYFy5EYqsLAOG9u/Mfd65F5lLK/PYuoFVkZKaBlbC/VyfjujVuRX3Zy2+EYO64+ulm7h9lVKJOoFVuVYuzk4O59XegVaBngX/KDo+16kwliDZXkdi9bME82ZfcqqhN8dkMABREVmIefgUdHyrqDjIq1RiKVT67FLJjXA51Rl0+OXSRpwOP5nPqFqpwgtNhqNzza6FdkglN4WWTNIBVHIjqdwOOXN1VsNJpUAibVw4pF9hBklRcpN2PhlPziTZDKvCK1Wh0Dp2x4Dx0LL1uzFxTH9+UAjbxLx6M5QD/sVrdxEdmwjf8l5o1rA2VColB/574Q9x6tw1dOnQDBXL+2DH/hN4Z/IwDtUGgxGh4VHwr1whl102CeZr3qIN/IWYM6eOxdFTl7B2y35MGNMfnds2A7/DpFBg2fpdqF+7Ojq1bZpv7jnjdXHWYt/hv3Djzn28PWkYvvr+F0THJuCjt8bgwcMYLF27A62bN+Al1Ox5SNYsvfiUzZdBNWubdvyBxvUD+SEpC1dsw6jBPeBb3jv77+wiY878dZg8biCSU9K4VpbGMn2On73CLwDee2Mk1CoVvlyy0Wp8OX3ZTLiAHQr9plhrL5wSMCa7Tdm6wJDTKTd2T4o62lSAgN6mRHZ3OBV+EuvOr7Han0H9B10/hr+nv902WUcC+kLJJXlnAnrJJS+0QwL6Qkvm8AApgD7lj3hk3HpiM0af5ytBXcmx3VzGRJ98tRIMjhl0JiSl4JmmdTH6+Z4Ivf8QPl7u2Lr7KAL8K8Pbyw17Dp3BpLEDsWbzPgS1bcL7sJ36vKCc166XpxtGDemBTTsOw9PdFZV8y+Fu2AO+I92wXo3s8fosA+b/tAWD+wWhTs2qBQI9uwB5kpaB/yxYh3cmD8eG4IOoWsUXGZl6JCQmw7eCD98ptwX085dt5RcgrEU+isX0iS9YBXqVSoVvf9yMd18bwftaGtuwbgAWrwpGz07PIHjvMQ78tQL8+EWBtfg+/24tUtPSUaNaZYwc1B0+Xh428y5EB6svlnrjo/mYPeMlLkTOxq742K2P7+dN51c6xdlYLBNmfIXk1LRcYQzq0zG79EYu59AXp06lzTcBvXAZnX3wX4hOfVygwU41OmFU87GFckpAXyi5JO9MQC+55IV2SEBfaMkcHiAJ0B9NQMZ1y/XzOQMvKtB/v2YHhg/syqH+VkgEP9zkhQFdsWXnERw69jcvjRnYuwMvyWG7103q14S55CYn0LOSkr2Hz6Bj6yaYOXUMctpVKpVw0qix8pe9GNIvCLO/Xo1xw3rnuiBg5Yxrtuznz7uNe6G3xdLGvHcU2L/Zjjnbod/46yEM7tsJC5ZvRfNGtVHVz5fbzwn0bVs0xKdfr8Tl66F4ZWQ/Pp8zF65jYK8OXNKdv58E65N3h54Bd1LKE4SEReKz915BkwaB/G6GpbGVfH343Jn9M+dvIFOnx+ghPTjQW4uP9WH6/H35Nr8I+Nc7L0LrpHF4fdo7sETv0Ns7SWv9pDzlpqix0vh/FCCgF2Y16LIy8fauaTCh4CdXa/rUxHtdPiyUUwL6QskleWcCesklL7RDAvpCS+bwACmAPv1aKlL/TLAZo+/EqoBa2JIbVsd9+u/reP3FQThy6iIHY41axXfN27dqlA30apUap85dxYzXRvCSHAa5O/Yf57vcOUt52CRywnhEVAwCqlbC+u0H+Q5/UJumWL1pLypW8MFzvTtYfU4pL9Bfvh6C4H3H8fG0sfh+9a+8RCY9IxNeHu44f/V2PqAvaskNKytfvXkf37wOj4y2WHJz/8FjXoPfoE4An/PZizfB3n20fMNuq/GZk/woJh6LVwbzkiR2ASV2s1pDn/clTSwQ8474hDEDso+CFDtAMe0T0Iuprni2CeiF0TYzKxPvENALI2YJs5IP6E1A+hMd0pJ1cHJWw9VLC5UEJ3+UMNkkDZeAXjq5pQB6U6YR8RsfwphutDoxlyYecA9yvPLBWg19lw7N8cmXK8F2mxmgdu3QHB1bN8Zn366BX+UK/GHQtycP4yC+fvvvOHziAmpWq4JHMXG8bn3U4O65SnlyltyY6/XZpMw7/ezvS1b9isb1avK59ujU0moNPSsRctJocC88CoEBfvjgzdFg9vPWoTMmZRcifbu1wZdLNnC7rK+5vMbMqJbq4OMTk3k5TY+glnh1dH/8uO637Hr6A0fP4dylm+jfoz1u3r2fq4a+fp0AXtfPduTZce3srsO3P27h8/nz9KVsG8y3Ob7BfYPw6X9X8QeKQ+9H4aURfbMfpBV7RRd4yo2lkpbFc6flq68SO0ix7BPQi6WsuHYJ6IXTl0puhNOyJFnKCfQPqzbA6Z23Ef/on3IAtUaFtgPqoEqgT0maVqmKlYBeunRKAfRsNvoHGUjcGWNxYmpfJ7ByG4h0ii6rZ9frs3LtFFv6jAXHwJWdQMN2lYU6BYyVRmf+/2k1OZubqzN/eLc0Nlbv7+ys5XdCpGql4thKR8UioHdUueIdR0AvnP62HorVqJzwYdePUMXDr1BOqeSmUHJJ3tkM9OGbDmD3TS30GZaPRqzfrioata8KhdKxEgDJJ1aKHBLQS5dMqYCew7LOhLSLydCFZ8CQrIemijNcGrjBqUbpfkt8WQN66VbvP54I6C2ozn64rcF+cSSJfOZWgIBe2BWx4uwy/B15zqLRkc3GoHPNzoV2SEBfaMkkHWAG+oMz1+K6yfoJRiqNEj3HNoVHudING5KKb6czAno7hRKgm5RAL0C4ZIIUsKiAVaA3HwnJXjLAjiVa/s37qFm9MqbMXADz21hLuqa0Q18yM0hAL3zeQuNDcTT0MG7F3oKL2gX1K9ZH7zp94ONSziFnBPQOySbZIDPQ//LiEkT71S/Qb+NO1VG/TeGOLZVsIqXYEQG9dMkloJdOa/IkngIFvliqVg1//jDEh3N/wuSxz/Gjf+R0Dn1RZSGgL6qCxTOegL54dC+MVwL6wqglfd/CAH3VeuXRbkBd6YMs4x4J6KVbAAT00mlNnsRTwOY59GxXPifQy+kc+qLKQkBfVAWLZzwBffHoXhivBPSFUUv6voUB+nqt/dGkc3XpgyzjHgnopVsABPTSaU2exFOg0EBPO/TiJYMs26cAAb19OhVnLwL64lTftm8z0G+dsBRRvnUKHNB2QF1Uq1fetlHqIagCBPSCylmgMQJ66bQmT+IpYLXkhp0neuKvK/hy1mR8Pn8tL7mpWMEb496ai2EDu9I59OLlhCzbUICAXv5LhIBe3jkyA/31H3/DwTA3q8FWrOaFTsMaQkGH3EieUAJ66SQnoJdOa/IkngIFnnLDduPZQ7A5G51DL14yyLJ9ChDQ26dTcfYioC9O9W37znkO/b67Loi8E59vkEarRp+Xm8HZzcm2QeohuAIE9IJLatUgAb10WpMn8RSgYystaEvHVoq34ISwTEAvhIri2iCgF1ffolrP+6bYtORMPLgdj7jIZLh5O6NyoDfY7jy14lOAgF467QnopdOaPImngM0aenayTWlt9FBsycwsAb3880ZAL+8c5QV6eUdbNqMjoJcu7wT00mlNnsRTgICedujFW10iWSagF0lYAc0S0AsopgimCOhFEFVgkwT0AgtagDkpgT7zTiJiV1yGLjwZhmQd1JVc4damCsqPawSlq1q6SZOnUqdAgefQ9+rSCt06tCh1kzZPiHboS2ZqCejlnzcCennniIBe3vlh0RHQS5cjKYDemJaFuFVXkLgzBDDln5uqnDMqvtYc7l2qSjdx8lSqFLAK9Oy8+R9//g1fzJwEVxfnUjVpAvqSnU4Cevnnj4Be3jkioJd3fgjopc2PFECvu5eE+6/9XuDEXFtUhP8XnR2efFp6Jjb/dhjDn+sGVxctboc+wO2QCAzo1b5Am+xUw8CAKujavjnOXLiBZT/vQkxcIny8PTBsQFd0D2qJOf9/2uHD6DgooED5cp6YPnEYFq3cjn7d2qJL+2aIS0jG+3OWYsKo/qhSqTz+9d+VSExKRfPGtfHe6yPh4e6aLwYWr9luVpYRr41/Dh1bN8apc9ewcOV2uGifPpD/fP/OqBtYDQuXb8Ost8ehsm85nD5/HT+t24n//ut1eHq4YUPwQezYd5z3H9Q3CKOH9ETkw5h8cTg7azF/2RYcO30ZTRvWwgdvjuKxZWTqsG7rATz/bGeU9/Hkdr5euglXboRCpVSiTmBVPmcXZ/keEmC15IYdTxkRFW1xEVTzq4h1i2bCt7y3wwtPDgNph14OWSh8DAT0hddM6hEE9FIrXjh/BPSF06s4etMOvXSqSwH0CZtuInblVZuTqrVtEJTuGpv9LHVIS8/AsvW7MXFMf74RyzZmr94M5YB/8dpdRMcmwre8F5o1rA2VSsmB/174Qw7QXTo0Q8XyPtix/wTemTwMzlonGAxGhIZHwb9yhVx2mW/ma96iDfxI25lTx+LoqUtYu2U/Jozpj85tm0GpVEChUGDZ+l2oX7s6OrVtagHo/4mXwf03Szfjw7dG4/yV27xvzgoRNpfFq4LRPagFhvTrxPveDo3gm85nL95C1OPY7OPU2QWKX6UK6NO1db441CoV7tx7gPHD+vC5atRq9OjUEhu2H8Txs1e4PTPbLlyxDaMG9ygxrEun3Fj4qaBTbhz6XSLZIAJ6yaR22BEBvcPSSTKQgF4SmYvkhIC+SPIVarAUQP/467NI/v2+zbiqze8O5wblbPazBvSffLUSLs5aMHBNSErBM03rYvTzPRF6/yF8vNyxdfdRBPhXhreXG/YcOoNJYwdizeZ9CGrbhPdhO/V5S60ZvOe06+XphlFDemDTjsPwdHdFJd9yuBv2gO9kN6xXI3u8PsuA+T9tweB+QahTM38pUc4LkEcxCRyq33tjJE6eu4qNwYc4lLPGgLuctyfOXrqJh4/j0KpZPTyOScD122F4e9IwftEwqE8QzIe4MPjfffAU3n9zFB+fM46T566hun9FHiO7cLhxJxxjnu/J++UF+M+/W4vUtHTUqFYZIwd1h4+Xh0N5kWoQAT0BvVRrTTA/BPSCSSmaIQJ60aQVxDABvSAyimqEgF5UeXMZlwLooxecR9KeUJuTKirQf79mB4YP7Mqh/lZIBK+0eGFAV2zZeQSHjv3NS2MG9u7AS3Ia1w9Ek/o1YS65yQn0Xy3ZiL2Hz6Bj6yaYOXUMctpVKpVw0qix8pe9GNIvCLO/Xo1xw3rnuiAwmUxYs2U/lAoFxr3Qm+/W523mCwW2mx/+IBrzP38LVSqWA3sHUnJKGjq0asSHuLk64174I363wb+KL1Zu3IOPp43lcVsD+h37j/M7B3njYDGbL1qu3LzHbbJdeEtAn6nT83n+ffk2gvcew7/eeRFaJ8funthMvAAdCgR6Jta3P27OdsOuxJZ/8372VZAA/ovVBJXcFKv8DjsnoHdYOskGEtBLJrVDjgjoHZJN0kEE9NLJLQXQJ+0KRfSi8wVPSgHU3jEECq3KoclbK7lh9d+n/76O118chCOnLnLw1qhVfNe8fatG2UCvVqlx6txVzHhtBC/JYTvdDIynT3zBYsmNubwnIioGAVUrYf32gxyWg9o0xepNe1Gxgg+e693BIsyzCeaM9/hfV3Hu0k2898YoHP/rMp9/3pIbBt+sPv7+g8eo5ueLOfPXcaA/eOxvuLu68AsV1nYeOIkMnQ6D+3bKF8exM5fx4GEMh/jf/zyHhMQUXpJkCejNSXgUE4/FK4P5BQK7EJJrswr0DObZFV3OWnmW3AkzvsLcmRNLxek3BPRyXZYFx0VAL/+8EdDLO0cE9PLOD4uOgF66HEkB9MYUPcIm7oMhIdPqxLwH14bv680dnrg1oO/SoTk++XIlKvn6cBju2qE5f/j0s2/XwK9yBf7g7NuTh3EQX7/9dxw+cQE1q1XBo5g4tG7eAKMGd7dacmOu12dBm3f62d+XrPoVjevV5HNhJTO2aujZHYXVm/bx/gHVKuGXX//IV3LDgH7EoO68D5urGehZvf9X3/+C9Iyn2jJb778xEucu38oXR/NGtfGfhT/zOwfxiSmY/e5LvJSGlfvsOngKjevX5GVIzOan/13FHwwOvR+Fl0b0Ra/OrRzOjRQDC3wo9oMpo/KBO7sVsm7LASyeO63En35DQC/FEhPeBwG98JoKbZGAXmhFhbVHQC+snmJYI6AXQ1XLNqUAeg6h56MROfNPi8dWauv6oPp33QF1/tIUIZRgdeR6fVauHWZLnzFfrEyFPaTKdqMtlco4Ek9yahoyM3W5hrJSGiFPUWSQr1KpbJbFFGZ+KalpYCfjsDsacm+FfrEU26Wf/c1qfD9veol58tdaEgjo5b48LcdHQC//vBHQyztHBPTyzg+LjoBeuhxJBfRsRuw8+oQtt/Dk3CPoo1Lh0rgCvPoFwq1dFekmXAyepAD6YpiWrFxaBHp2lTNl5gL+kEPep50J6GWVvzIZDAG9/NNOQC/vHBHQyzs/BPTS5kdKoJd2ZuStLClgtYbeWmkNq5EKCYvE3I8mlnidaIe+ZKaQgF7+eSOgl3eOCOjlnR8CemnzQ0Avrd7kTRwFHHqxVM5QSvJLpgjoxVlUYlsloBdb4aLbJ6AvuoZiWiCgF1NdYWxTyY0wOtpjhYDeHpWoj9wVoHPoLWSIXiwl72VLQC/v/LDoCOjlnSMCennnh3bopc0PAb20epM3cRQgoCegF2dliWiVgF5EcQUyTUAvkJAimSGgF0lYAc3SDr2AYtowRUAvndbkSTwFCOgJ6MVbXSJZJqAXSVgBzRLQCyimCKYI6EUQVWCTBPQCC1qAOQJ66bQmT+IpQEBPQC/e6hLJMgG9SMIKaJaAXkAxRTBFQC+CqAKbJKAXWFACeukEJU/FogABPQF9sSy8ojgloC+KetKMJaCXRmdHvRDQO6qcdOMI6KXTmnbopdOaPImnAAE9Ab14q0skywT0IgkroFkCegHFFMEUAb0IogpskoBeYEFph146QclTsShgFehj4hIx7q25iIiKzhdYST6qMudk6NjKYllzRXZKQF9kCUU3QEAvusRFckBAXyT5JBlMQC+JzNyJlDv08ffS8feGSCRFZiAzNQvuFZzg38ILzV6oAo2LUrpJk6dSp4BVoJ85bxmfbGl4gZS1rBHQl8z1TEAv/7wR0Ms7RwT08s4Pi46AXrocSQH0+nQjLvwShVu/xwCm/HNz8Vaj9fhqCGjvLd3EyVOpUsDqi6Xe+Gg+Zs94CY3q1ShVE6Yd+pKfTgJ6+eewNAL9Q1UydmqvI1KVhBRlBryMLqidVQEDMxrC0+Qs/6TkiJCAXv7pIqCXLkdSAH1ieAZ2fnCjwElVbuyBXh/XdnjiaemZ+PTrlWjeqA7GPN8Tt0Mf4HZIBAb0ao+wiEf4YvEGPIqOh7NWg5GDe6B/z/bYEHwQew6ehk6vR83qfvjordHwLZ//ouLE2av4Yc0OJKemwcPNBTOnjoVGo8bC5dsw6+1xqOxbDqfPX8dP63biv/96nftlrX2rRvzPVb/sxXN9OmLrrqM4cPQsendpjcnjBuaa66lz17Bw5XaolErUr10d0yYMhYe7q8XYB/cNQtSjWHz85XK8M2k4mjQIzLaVkpqGL5dsxPhhfVA3sKrDepa0gQT0FjJGL5aS9zImoJd3flh0pQnodTBgr/MN/KG9C6OFrTUXk4ZDfZAuEAr5p4ZHSEAv/0QR0EuXIymA/uqOx3yH3lYbsbwpnNxUtrpZ/D4tPQOLVgaDAe2bLw9BfEIyrt4MxbM92mHOgnWY+upQ+FeuAIPBiNDwKNwJjUTU41i8OupZKBQKxMYnQalUopy3Ry77CUkpHJA/eHMUfLw8kJGp4zbYRcLiVcHoHtQCQ/p1wjdLN+N2aAS+mDkJV2/d4za6dWjB/1y4YhtGDe4BN1dnLFu/GxPH9IerS+6NkMMnL2SP2RB8CJ7urujSvpnF2OvUrIrte/5EeGQ0XF20mDT26cWBPsuAJauC8TgmAWOH9irVm9J5F0GBJTe1avjzRJfWRiU3JTOzBPTyz1tpAvpkRQbmeBxEukJvVfgqBk98kNodqhKC9AT08v8ZIqCXLkdSAP2JH+4j9M94m5Pq91ldVKjjZrOfpQ4M6BksN6wbgJCwKAS1bYJrN++hcf1A7Nh/nO+q52xzF/6MQX2CbEJvpk6PWV8uh39lX/QIaom6tapB66TBtVthOHvpJh4+jkOrZvU4RF+/HYa3Jw3jQL8x+BD8KlXgLiMeRnPQtwfoO7dthqVrd6B18wa8v6XY2UXFdz9twchB3bFi4x7MeG04vD3dsXX3UTSoE5A979JcZWI30LNE/fjzbzwBea+iHFppMhxEQC/DpNgREgG9HSIVc5fSBPR/ax5gtetZm4rOSu2JSobcO1s2BxVTBwL6YhK+EG4J6AshVhG7SgH0p5dH4M6hWJuRCgH0Lw7rgx/W7kCzhrWRkvqkyEBv3vm+FRKOvy7cwIEjZzFtwgvw9HDjdwD8q/hi5cY9+HjaWA7XZqBPTklDhxwlNy+P7JcL6C9dC8Hsb1ajcsVy+Hb2m7h8I4RfBLDSoW4dW+CVkf1w/fZ9i0B/594DBO85BmZz46+H0KJxXb5Tz+48dG3fHDv2n0CjujXQpkUDqFRl42FjqyU31k64YYmlU25s/kxSBxEVIKAXUVyBTJcmoN/jfBN7tQXXvjLZXkprjWf0JaNek4BeoIUuohkCehHFzWNaCqC/fTAWZ1Y8rSu32hTAqFXNoNY6BqDmHXpWzsLq5xcs38pr1Xt1boV5i9bzHfryPp7cfVLyE+z54zSUCgVGDOrOP9Pp9DAYTXBxdsoVIitj0ev12Zu7v/95DlFsV75pPQ70g/oG4f6Dx6jm54s589dlAz0z4kjJDbPLQP+l4X1RpVJ5i7Fv2XWElwlVreLLS4XCIx/j5RH9cOl6CI/95LmrqBtYDUP7d4azNvd8pFtZ0nqic+gt6E019NIuwsJ6I6AvrGLS9y9NQH9Aews7na/bFJGA3qZEZatDVhaQlACkPwE8fQBPr0LNn4C+UHIVqbMUQK97YsCOGdeRkZRlNdb6fX3R+kXHNwVyAr2TkwaffbsGAVUr89LpY2cuY/HKYF4uExOXwMtZhg/syuvT2Y44e9A1PjEZM6eOQ3X/irliTExKxYx/fw9XV2fej9Xofzx9HOLin9bomy8ImH9bQP/H8fPYvPMI920eZ3aWs4aePfA6b/EGfDRlNELuR+WKvUmDWvz5gKkThsLLw43X9LPyIVZHzwCftU07/uB3Jqjkpkg/GiVnMJXclJxc5YyUgF7+eStNQH9J8xDLXU/bFJ1KbmxKVGY6KK5dgmrHRkCXmT1nk391GEa+Cng83SG11QjobSkk3PdSAD2L9uGVFBycd9fisZXlA13R97O6UKrEfbSewTvbgWe72+bGduBhMvFTaxiUP0nLyCWuVuvEH1B9ulOfxUtbiqNZir044pCrT9qht5AZ2qGX63J9GhcBve38qJL0cIpMhSpZD6OLGroqLsiq6GJ7oEA9ShPQZ8GIr9wPgx1baa010VfBpLR2AqknvhkquRFJY6MRyv2/QnnmmGUH7p4wvDAOphq2jyYkoBcpR1b+z5fKGzuP/tqux4i6mIyUx5moWN8ddbqVR9VnCncHR6x4CwJ6sXySXWEUsAr07KHYCTO+4meO5m1UQy+M+GTFMQUI6AvWzflmErz3RkChM/7TUQGkdKmM1DYVIcVBLKUJ6JmIUapkfON+BOwIy7ytnNEVH6R2g6up5NRpEtA79rvH1ihFbDSUP34Nhd76iUimwLowjH/dlil6sZRNhYTrINUOvXARkyVSIL8CFoGeXaFNmbkAHds0QbuWDXOddsPeINurS6vsBx1KsqhUclMys0dAbz1vzneS4bM9zGqHlI6VkBpUSfTElzagZ4LpYcRRbQguaaLwUJmMGoZyaK2vija6ACmukQTNGQG9oHJmG1MePwjlwd02jWfN+LfN0hvaobcpo2AdCOgFk5IMFaMCNl8sxWJjTxt/P286f3sYe2hh3ZYDWDx3Wok/zpKAvhhXXhFcE9BbFk+RaUSlH25AkZl/FznniJiX6iCrkrjlN6UR6IuwZGU3lIBenJQot66B8upFm8YNYybBVKdBgf0I6G3KKFgHAnrBpCRDxaiATaCvWMEbH8z5EV/OmsyBnpXi5AT8Yoy9yK4J6IssYbEYIKC3LLsmKg0V1t21mZOk3lWR1qKczX5F6UBAXxT1xB9LQC+Oxsrfd0F54pBN44YpH8JUoeA7ZQT0NmUUrAMBvWBSkqFiVMBmyQ077oiV2ZjfGsteGnDiryu0Q1+MSSvrrgnoLa8A1wvx8DrwwObySGtWDkl9HT8azaaD/3/NNwG9PSo51scE4GpYFh7EGJGcZoSrVgFfLyWeqaOBRm2fTQJ6+3QqbC/FtYtQbVlT8DCNFlkz5wE5ThmxNICAvrDqO96fgN5x7WikfBSw65SbmLhEmF80xY4uWv7N+6XibE/aoZfPQixMJAT0tENfmPVSmvrGJhux92wmImNzPPD8vwl6uCrQv40WNSqpbE6ZgN6mRI51MJmgWj4fishwq+ONg0bB2KKNTfsE9DYlEqwDAb1gUpKhYlTALqAvxvhEdU1AL6q8ohknoLcsLa+hX3IdCn1+2Ms5IvalOtBTDb1o61NMw7vOZPLdeWuN7dSP7+ECjabgKAjoRcxSQhzUP34LZOQ/Ic7YrDWMQ0bb5ZyA3i6ZBOlEQC+IjGSkmBUgoLeQADqHvphXpQ33BPTWBXI7GwvPP6KsdpCi3IY5p5Ib4X+GUtONWPxbuk3Dz7X37oZSAAAgAElEQVTTomFAwbU3BPQ2ZSxyB0XYXShuX4ciPhbGgECgYTOYvHzstktAb7dURe5IQF9kCcmADBSwCvTmoyvPXLjB3xDGymxqVq+cfZwlq60v6Y126EtmBgnoC8ibEaiw7g40j/KDn9FdjegJ9WHSKkVPPAG98BLfjzZg4+Hcb3C05KVjYyd0alTwFj0BvfD5EdoiAb3Qilq3R0AvndbkSTwFrAK9+UHYUYO748O5P2Hy2Od43TwdWyleMsiyfQoQ0NvWSZWgg1PkE6hS9DC6qqGv5Ax9ZVfbAwXqQUAvkJA5zJwPycKBc5k2Ddevrsbg9gW/mp2A3qaMxd6BgF66FEgJ9NcSMvDfC9EISdYhMdMAfzcNuvq7Y0rjCnDXiL/ZIp2q5ElqBWweW8l25XMCPR1bKXWKyF9eBQjo5b8mCOiFz1FcshHL9touuendSouWtajkRvgMSGuRgF46vaUA+lS9Ed9eisGGOwlgJ1Xlbb7Oasx8phKere4h3cTJU6lSoNBAX9J26Nkxm9/+uJknrW2LBrmO26SSm5K5lgno5Z83Anrhc2QyAesOpiMq3vpDz84aBV7p6wJPV0WuAP7z71n4ftG32Z/1B7ALADtr5WyOnrNm/wevT3lb+ODJYqEVIKAvtGQOD5AC6G8lZuK5vfcKjLF9JVes7l7d4XmkpWfi069XonmjOhjzfE/cDn2A2yERGNCrPcIiHuGLxRvwKDoezloNRg7ugf4922ND8EHsOXgaOr0eNav74aO3RvN3DuVtXy/dhCs3QqFQKNCqaT28NKIv3N1ccOrcNWzZdRRzPniFv2x0+54/cerva5g+cRj++/0viHociw+njEbLJnWzTbI4N/92GMOf6wZXF212nKwKZM6CdbzfrGnjeKwLV26HSqlE/drVMW3CUKhUKotj2RxZY8z34GEMPp42FmqVymp8/373Few6eAq/7T8BJ40a5ct5YubUsfDxKrkXVFZLbsznzbMXSn0+fy0vuWEvmWLHVw4b2BUloYaeXXx8uXgj1i2ayRcoKyNibe5HE/mfBPQO/94o1oEE9MUqv13OCejtkqnQnRJTTVixPw16KwfdDA1yRh3//MdWEtAXWupiH0BAL10KpAD6H6/H8R16W+3sC3Xh6WDpDXv2cdHKYKSkpuHNl4cgPiEZV2+G4tke7TgoT311KPwrV4DBYERoeBTuhEZy4GY8x0A9Nj4JSqUS5bzzQ+3CFdswanAPzlIXrt7B1l1H8cnb4zm8bww+hPHD+qBZw1qYt2g90jMy8Z8PJ+Be+CPuf8Sg7rmmzeJctn43Jo7pzy8CWOWHuR/jNta6dWjBS7zNf98QfIg/z9k9qIXVsQlJKViy6ldkZRnw6uhnUc2vIrdhLb6dB06icf3AUnEMO9OpwFNumBBTZi7IlYjFc6dxoUtCy/lCLBZvXsAnoC8JWcwfIwG9/PNGQC9ejtIzTfjjkg63HxiQqTdBpQSq+arQs6UTKnharsEloBcvH2JZJqAXS9n8dqUA+g9OP8Sv95JsTmpTrwA0r+Bis5+lDmZQblg3ACFhUQhq2wTXbt7j0Lpj/3G+A52zzV34Mwb1CbILaHMCvclkAhs7tH8XPIyOQ0xcEkLCIrkddrFw9uJNzJo+rkCg/+SrlXBx1vJddAbizzStyy8KLAF957bNsHTtDrRu3gCN69eAtbHHzlzBw8excHdzxZP0dAwb0JXbsxYfuyj5+/JtVPevxDeqq/tXdEh3uQwqtcdWmk/p6dimSfbdBHYV+M7sJfh29pt84RHQy2UZFi4OAvrC6VUcvQnoxVfdaAJiEg0o56Gy+YZYAnrx8yG0BwJ6oRW1bk8KoP/X2UfYdDfR5qSEAPoXh/XBD2t3oFnD2khJfSI40LNJmC8GouMS+JzY7v6fpy9jxuTh+OnnnfmAftOOP7Dyl71o2jAQH7w5Gqs378PwgV051N8KiUBEVDS/A5AX6NnuOivR6daxBV4Z2Y/v/n+/Zke+sSMHdceildvRrmUjuLs5g+3os4sKdgfBWnwatRpqtQqh4Q+xeOV2fPrOS/D2creZI7l2KPVAP25Y7+w7CnmBvriT0q1bcUdA/kkBUqAsKBASEoKIiIjsqfbHaezCR2iDH3AW9bM/r1WrFqpVq1YWJKE5FrcCJsCQooMp04Cdnz2ApoobXFtVhsLBcpPino4t/xv/H+Znn31UYDf25MuF4fXgosr9DIwt2+bvc5aysPr5Bcu3oneX1ujVuRUvhWE79OV9PHn3pOQn2PPHaSgViuySGJ1OD4PRBBdnp3wuc+7QRz6KBfs3q3M/d/kW79uycV0kpaRyQP/upy02d+gLU3LDavZnf7MaLw3vi1o1/CyW3LDN20UrtqN7UEsez5GTFzF2aC+YLzgsxcfKfVh7kpaB/yxYh3cmD0eFcl72yi27frmAPiYukdfIsyulghqrSzLXpctuRv8LyJ4d+uKOnYC+uDNA/kmBsqEAAX3ZyHNJmaVJb4QuLAnGtKcPg6xrfZT/6VzPB36zO0JV7ilolaaWpDPi2d0hiM0wWJ3W+Ho++LhlJYennRPonZw0+OzbNQioWplXKRw7cxmLVwajbq1qiIlL4OUrbIec1dazHXAPNxfEJyZj5tRxFktPPv9uLVLT0hH1KA5urlp88vaLvB4/5446C5xxJAP6qROG4vvVv/Kym5lTx6BBnYDseVmroQ9q2xRfLtnA+7Fd/NuhTzchWJl31KNYzFu8AdMnPI9dB8/kq7/3cHeDCSb079GOj2HlN9dvh6FurarZNnLGx+CdxcnKh9jDt107NMe4F3rzZwlKarMI9EnJqfxFUqwspSQ3qqEvydmzHjuV3Mg/r1RyI68cUcmNvPJhTzSltuQm04C094/CFGP5CFaFrwtcv+wCOOd/uNse3RzpI0XJDYvr5KMneOVwhMVjK5uUc8bGXjUg9g0KBu9sBz4nuOqzDIDJBI1GDQbbbMc6Z9NqnfgDqaWxsfmyk3O0TgW/jK8kzN1iyY35qMe8xzyWhAnljJFOuSlpGbMvXgJ6+3Qqzl4E9MWpfn7fBPTyyoc90ZRWoNetvQb9vrACJdD0rQGn8Y3skUmQPlIBPQuWnUe/4kY8jj1Mxf0UHVpVdMWwWt7o7i+P2u2yBvSCLCCZGLFaQ28uWTlz4QavKyoJx1Ra0pTOoZfJShMwDAJ6AcUUyRQBvUjCOmiWgN5B4YpxWGkF+rQZR2B6+KRAZRUBnnCd10ky9aUEeskmRY7KnAI2H4plD5JOmPEVF6Y0lOHkzDCdclMy1zsBvfzzRkAvrxzN/ewTLFn4TXZQ1l4s9cm/5+K1N6fLK3iZRWOCEVnKBOgVSdCYvKAxlmMnQAseZakF+hf3gtXQF9QULmq4rugjuKbWDBLQSyY1ORJRAZtAn/Ms+pJ0Br09mhHQ26OS/PoQ0MsvJ3kjIqCXd458jv4Ol6GD8OTPk0hq2DxfsPonBoQfjEdajA6GLBNcK2jg18Eb7n5aeU9M5OiS1X/hsfMWZCn+OU9cbfRBlYzRcDc0E9R7qQV6e3bo2Yk333QVVM+CjBHQSyY1ORJRAatAn/PEm0F9Oma/XVXEWCQ3TUAvueSCOCSgF0RGUY0Q0Isqb5GNFwT0SWEZuLriATIScr+OVqlRoOF4P1RqIY9Xo5vY2deZBkSnZ8HbSYWKzmpRHyiM0xxEtPNmq9r7ZbwIL33HIufGbKC0An1Zr6EXbIGQIVIgjwIWgZ6dDrNj/wmU9IdibWWbgN6WQvL8noBennnJGRUBvbxzZA3odclZOPNFGNiflppSpUDLtwPgVaN4jxV8nGHAzPPRCE3VZYfpoVHis2a+aC7CkYc6ZSxC3WaBldtYawo4oXbqXKhNT8/5zttiYs4jIvwQ0tMeQ61xg6dHAOrUGw212rKWpRXokWFA2gfWT7lR+rvDZU4QoC19p9zI+7cCRVfSFSjVx1baSg4BvS2F5Pk9Ab0880JAL/+8mCO0BvR3tj1G+OGnb3601rwCXdDqnX/OlJZ61qdj0jH3aixSLNRhKxXAhDreGFnDS9CqdlZqE+my3OZUq6VPhXtW4zz9TLh5Yx2uX10Gkyn3BYGXVy207fAfeHjkf5lXqQV6AKaEDOh+voGsU1G5tFK394PT2AZQ+Eh7wUglNzaXNnUoAQqU2hdL2aM9Ab09KsmvDwG9/HKSNyLaoZd3jqwB/bmvw8BKbgpqKq0SXb6ui+J4/4rRBLx55iFuJv+zM583Vg+1Ess7+KGigOeYx2h3ItZpp82kVswYjvL6nrn6RYQfxF+nP7U6tly5hujc/XuolLnPwS7NQG8WwxiXAePVGP5PZWNfKMtLC/LmOAjobS5t6lACFLD5UGwJmIPDIRLQOyxdsQ4koC9W+e1yTkBvl0zF1ska0B999zayMgo+gYQF3e6TmnCrJP0DsuFP9HjxRO5dXUsiTm9YHoOqCneud4LTUTzSrreZL//0CfDMapOr357fBiE9I7bAsS2feR81aw3K1acsAL1NQSXqQEAvkdDkRlQFCOgtyMt+uK3BvqjZION2KUBAb5dMxdqJgL5Y5bfp3BrQX1kViei/Uwoc7+SpRqe5tW36EKPDHw+f4PMrBcMx8zuwqjveaVhesBDSVWEIc51r015g6ufQmipl98tIj8XunblB3ZKRwMBBaNHqfQJ6mwqL04GAXhxdyaq0ChDQE9BLu+IE8EZAL4CIIpsgoBdZ4CKatwb0j/9OxtVVBe+AV+vqg7ov/AOtRQylUMOvJWZiyl+PbI5hdfRjanrZ7FeYDlHOq5CkOWV1SDl9D1TKGJHr+0cPT+HEsXdtumFlN916LiOgt6mUOB2kBPqYmBs4cXoB4hPuISMjCR4eVVAjoBPatpoEJyfh7iqJoxRZlbMCBPQE9HJenxZjI6CXf8oI6OWdo4KOrTy/IBwJd9IsTkDrpUb7TwKhclYWywTTs0wYdCQcNt5LhP8+UxGtyrsIGqMRmQhz+wKZysh8dl0MAQhI+xAK5D6ZJSsrHTu292KPgRYYS736Y9G46esE9IJmzH5jUgC9TpeKU2eW4PJVdvRp/vXg6loBnTu+izq1e9sfOPUkBXIoQEBPQF/ifiAI6OWfMgJ6eeeoIKBnNfShO2MQ8WdCLu4o38gN9UdWgbOPulgnt+V+Mr6/Zf0kno6+LpjToqJoMWYow5GiuQj2JwN5j6xnoDX6WfV35uTHePDgiNXvlUo1uvX4Cd4+9QjoRctawYalAPq4uDvYsDn3HZy8UVWr2gaDBy51WIW09Ex8+vVKNG9UB2Oe74nboQ9wOyQCA3q1R1jEI3yxeAMeRcfDWavByME90L9ne2wIPog9B09Dp9ejZnU/fPTWaPiW984Xw9dLN+HKjVBkZRkwtH9nDOnXCXfuRWLuwp+h/N8T8u2eaYixQ3vjs+/WoF+3tujSvhniEpLx/pylmDCqP9j3+w7/hdWb9iE9IxOd2jXFmy8Nhkql4nYuXLmDcj6e+ODNUWhQJwB/HD+PxauC4ebqgplTx/DPTCYTDhw9B28vd360Omubdx7Bb/tPwEmjRvlynpg5dSx8vOTxvgyHk+nAQAJ6AnoHlk3xDiGgL1797fFOQG+PSsXXx9abYllkhkwj0uP0MGaZ4FJBA42rdOeC21JmxrnHOB+f/zQeH60Kazr6gZ10I5em0yXj931jkZERZzGkps2nok7d/KBHD8VKl0EpgP7v86tw8swim5Oa9MpRaLWOwWhaegYWrQxGSmoa3nx5COITknH1Ziie7dEOcxasw9RXh8K/cgUYDEaEhkfhTmgkoh7H4tVRz0KhUCA2PglKpRLlvPP7X7hiG0YN7gEPd1fMmb8Ok8cNRHJKGrc/YlD37HmxGOYt2sBPwWJgffTUJazdsh8TxvTnkL374Cm8/+ZoaNQq7P3jDMIjozFp7AAYjSaoVEpcuXkPew6d5p+xubz72nDExCXxC4/33hiJU+euYd3WAxg1pAe6dWjB/W7a8Qca1w9Eo3o1bOpbmjsQ0BPQl7j1TUAv/5QR0Ms7R/YAvbxnAESmZeFYdBquJWYgwM0J7Xxd0Nhb+pN37NHJkJWBq1d/QsT9/cjMTIRCoYRPuQZo2uwtlK/QxKIJAnp7lBWmjxRA//sf/8LNW7tsBjzs+dWoXKmpzX6WOjCYXrZ+NxrWDUBIWBSC2jbBtZv3OOzu2H+cA3bOxnbFB/UJsguEzUDPdtO//XEz3n1tBCIfxWL+sq38IoG1Jg0C0a97Gx6Dp7srKvmWw92wB1AplWhYrwZC7z9EYECVbBCPiUvEguXbeFyuLk9/dhnkswuSpg1rIXjvMXz01hgkpTzB0rW/4a1XhsDVxRmHT17gfc1Azy4Y/r58G9X9K2HYwK6o7i/eHTqHEiPRIAJ6AnqJlppwbgjohdNSLEsE9GIpK4zd0gD0wightRUTkpJC4epaCRpNwQ9AEtBLlxspgP7w0f/g6vVtNiclBNC/OKwPfli7A80a1kZK6hNBgP7z79ZysA4Ji8Rn773C4f3arTCcuXAdA3t14PPSap2gVik50A/pF4TZX6/GuGG9s0HeEtB/99MWzJo+joM6K+n5efvvmDVtHN+5N1+EPEnLwE8/7+R3BSwBvV6fBbVahdDwh1i8cjs+feclXpJT1hoBPQF9iVvzBPTyTxkBvbxzREAv7/yw6AjopcuRFEB/5dpWHPnT1tGnCrw+8QTUasdesGXeoZ84pj+vn1+wfCt6d2mNXp1bYd6i9XwnvLyPJxc2KfkJ9vxxmte/m0tmdDo9DEYTXJyd8olv3qFnNfirN+/D7Bkvcei2VHLDgJ7FEBEVg4CqlbB++0G+M69WqXHlZiheH/8cL/Fhu+rHzlzG9Ikv8FKaI6cuYtqEoRzaWfnPD2t28DKbxORU/LhuJz6cMhpaJ02+HXpzsAz8/7NgHd6ZPBwVygl7ypV0q9FxTwT0BPSOr55iGlmWgD7GlIgN+B1heIQM6OCPCuiHdnhGkfsBumJKhVW3BPRyy0jueAjo5Z0fAnpp8yMF0GdmJuPnjUORlm75WQo242ZNRqFz0HsOTz4n0Ds5afDZt2sQULUyr5Fn4Lx4ZTDq1qqGmLgEtG7eAMMHduW19exhWg83F8QnJmPm1HEWS1bMQM8emGUPpZ67dBP9e7TnD61aKrlhQM/AnLUVG/dwoA9q0xRLVgXz+n13VxdeWvPx9HG8JOetjxegVg1/aNRqfhHAauS37DyCC1fv8n6jn39aM8/mwWroWRv3Qm8807Qeh3j2sCy72OjaoTn/nF0wlLVGQE9AX+LWfFkAeiOM2Gg6hEWmrUhB/iMEn1d0wQzFCHjDsYenxE46Ab3YChfNPgF90fSTYjTt0Euh8lMfUgA98xPx4Ax+3fmGxWMrK/o2BCu3YaceidkYvLMd+JzAq88yACYTNBo12EUB2+nO2VgpDauJF6plsjsBBkM28Bdkl8XL4mIP0RbcL4OflsN28MtqI6AnoC9xa78sAP1dPMAI46d8V95a+0AxBuMVfWWZPwJ6WaYlOygCennnh0VHQC9djqQCejYjdh79+YvrcD/8BJKSI+BXpSUaNRiMmjW6SDfhAjxJAfSymGgpDIKAnoC+xC3rsgD035o2YYWp4BMR/FAevyvnyzJ/BPSyTAsBvbzTkis6AnrpkiUl0Es3K/JU1hQgoCegL3FrviwA/STTVzhhulJgbpRQ4rTiR7gpHHuASszEE9CLqW7RbdMOfdE1FNsCAb3YCv9jn4BeOq3Jk3gKENAT0Iu3ukSyXBaAvpdxOqJg/eEps7TbFJ+jvkJ+L9MgoBdp8QtkloBeICFFNENAL6K4eUwT0EunNXkSTwECegJ68VaXSJbLAtD/27gSm3G4QAXLwwt/KheLpHLRzBLQF00/sUcT0IutcNHtE9AXXUN7LRDQ26sU9ZOzAgT0BPRyXp8WYysLQH8G1/GKcV6BuRmGrpitfFWW+SOgl2VasoMqq0CviNFDEZ4JZZoRJo0C8FDBUM8FUMvviDsCeul+hgjopdOaPImnAAE9Ab14q0sky2UB6Jl0X5h+xjrTfosq1oI/Nin/DRfI81X3BPQiLX6BzJZFoFfeSIPqahpgyi2iyUuNrA4eHO7l1AjopcsGAb10WpMn8RQgoCegF291iWS5rAA9k+80ruFr0y+4YQrjalZQeGMkuuMVxQBoId/zdgnoRVr8Apkta0CvjNRBdSLZqnoc6nt5A0qBBBbADAG9ACLaaYKA3k6hqJusFSCgJ6CX9QK1FFxZAnrz/OORjDgkow6qloh8EdDLO01lDeg1v8UDGcYCk2Jo6Q5jbfmcGEVAL93PEAG9dFqTJ/EUIKAnoBdvdYlkuSwCvUhSimaWgF40aQUxXKaAPt0Izc54m7oZazjD0MbdZj+pOhDQS6W0dG+KlW5G5KksKkBAT0Bf4tY9Ab38U0ZAL+8clSWgVzzUQX3MermNOVOmcmpk9fSWTeII6KVLhZQ79E9MmQjXxSLdpEOWyQCtQg1vlRuqaspDJaeaL+nkJ08CKUBAT0Av0FKSzgwBvXRaO+qJgN5R5aQZV5aAHlkmaLbbfqeDsb4LDE3dpEmAHV4I6O0QSaAuUgC9AUZE6OPwOCvRYtQaqBDgVBHlVfK5SySQvGRGIgUI6AnoJVpqwrkhoBdOS7EsEdCLpawwdssU0ANQn0yG4oGuQPHYQ7EmH7UwAgtghYBeABHtNCEF0KeZdLiScb/AiDyVrmig9bczasvdwiIe4YvFG/AoOh7OWg1GDu6BhnVr4PPv1uLjaWNRv3Z17Pr9FOrWqoa4+CQsXLkdKqWSfz5twlB4uLviqyUb4V/FF2Oe74mMTB1mf7MaLRrVxohB3XM51WcZsCH4IPYcPA2dXo+a1f3w0VujsWbLfly5Ecrturo4Y8Zrw7Ftz5/8s6wsA4b274wh/Trhzr1IzF34M5SKp8fGPt+/Mwb0bI8VG/fgwcMYHq9apYIlPx9OGYXjf13Bph2HkZ6pg1+l8vjk7fHwr1whnzCnzl3jn7Vv1QjhkY+xfc8xjHuhD75Y/DOaN6rD53k79AFuh0Sge1BLfPbdGvTr1hZd2jdDXEIy3p+zFBNG9UdcYjJW/bIXLZvUwTuTR8DF2alIuRJ6MAE9Ab3Qa0p0ewT0oktcZAcE9EWWUFQDZQ3ooTNC/XsiFE8sPxhraOQKYyNXUTUvrHEC+sIq5nh/KYA+KisBEfpYm0G2cqnlcOlNSmoa5ixYh6mvDuVgazAYERoeBZ0uC1t3H4WXhxumvDIE23YdReP6gYiOS+DxdOvQAhuCD8HT3RUDerXHl0s2Ijo2AR+9NYaD9dK1O9C6eQO8OurZXPHvOXQGUY9j+ecKhQKx8UlQKpX4edsBjBrcA77l/ylhW7hiG/+MXTDMmb8Ok8cNRHJKGq7eDM11oZCQlIIlq37l4P/q6GdRza8iLPkJj4zGjv3HMWv6eGjUKqRn6GAyGfkFRN52+OQF/lG9wGpYsGIb3p40jM910cpgMM3efHkI4hOSeSwDe3fAvEUbwK4xZk4di6OnLmHtlv2YMKY/18k8j5xzs5lUiToQ0BPQS7TUhHNDQC+clmJZIqAXS1lh7JY5oAegyDJBeSMditAMKDKfgr2xvBrGZu4wVZDPzrw5wwT0wqx1e6xIAfQhuseINdh+lqORthrclY6dtnTtVhiHXAaiORv7/PKNEETHJiKoTRPcvfcgF9B3btssG9rbtKjPobVqFV9kZOqRkJgM3wo+HHzzAj3bXR/UJwiN6tXI5Y/dDUhNS4eLVgutVoPxw/pg2+6jHOhVKhW+/XEz3n1tBCIfxWL7nj/RpkUDuLpo0bp5fZy9eAsPH8fC3c0VT9LTMWxAV76Ln9dPYlIqpn+6GF07NEdQ6yYIqFaZg72lxoD+4eN4XLp+F9NeHQq/yhWQlp6BZet3o2HdAISERSGobRNcu3mPAz37nAF/Jd9yuBv2gN9paFivBgG9PT9MxdUnKi7domv2w23tu+KKlfz+owABvfxXAwG9vHNUFoE+Z0YUSVkwuSgBJxkdPJ9nyRDQS/czJAXQ39NHIzoryeakxAJ6tvvMSk7Wbz/Id71bNK7Dd+g3Bh9CWnomunVsgVdG9uM77QzoB/fthAXLt6J5o9qo6ueL0PsPMbhvEN6ZvYSX88ye8RIYKFsC+vnLtqJvtzYo7+PJ7Xl5uuOLReuRlPIEIWGR+Oy9V9CkQSDYhcaZC9cxsFcHvrPP7iB8v+ZXtGvZCO5uzvyuwazp4/gFgCU/mTo9Ll0LwclzV3Hmwg3M+eBV1Arwy6cxi3Ppmt/g4+2O994YxfuYgf7FYX3ww9odaNawNlJSn2QD/ZB+QZj99WqMG9abzz0woAoBvc3VW4wdCOiLUfwiuCagL4J4Eg0loJdIaAfdlHWgd1A2SYcR0EsntxRAz2CeQb2t1tqlNpR4WlNe2MZKXuYtWs936BlMs5aU/ISXzZhLW5at34XzV+5gystDsktuWjWtx+vkXxrel4O2uawkPSMTXh7uOH/1NofavDv0G389xOvfzbX1Op0eBqMJy9bvtFpywy4GVm/exy8IWNlMzpIbFueiFdt5HTtrR05exNihvfjdhbx+UtMy4e6qhZPT0xcssnmx5wA6tW1qEejZh6zkZt7iDfhoymh4e7nznfiJY/rz+nl28dK7S2sM6tMx+/OIqBgEVK3EL4II6Au7GiXuT0AvseACuSOgF0hIEc0Q0IsorgCmCegFEFFkEwT0Igucw7wUQM9OubmUHgY9DFYnVlntjQCNb5EmfuzMZSxeGcwfeo2JS+C17+1aNswG54fR8XyXfda0cblq6KMexWbDLqu3z1kDz3a4LQG9uWaf7fB7uLkgPjEZM6eOw5rN+6yW3LDa8wNHz+HcpZvo36M9bt69n31BwGrlTTChf492XINjZ67g+u0wjB7Sgz8bkNPP+GF9MX/ZFsyFYW0AAB8vSURBVFT3r8T7ujhr8f4bI3mNft5mrqFnNfDswVx2QcH6/rLjMAd6dlHw2bdrEFC1MkYN7p4N9OZ6fPaQLgN69oAuu/vQI6glXh3dH9r/XUwUKWECDqYaegtiUsmNgCtMBFME9CKIKrBJAnqBBRXYHAG9wIKKYI6AXgRRrZiUAuiZ6yRjGm5mRlqMwk2pRSNtdQf35vObZPDLTmFhJS9iN3YKDUwmaDTiPotiyQ8ru2GNwTX7Pik5FSaTKXvKrJTH29MdKpV8y+uEyg8BPQG9UGtJMjsE9JJJ7bAjAnqHpZNkoCRAb0iHR8hPUKXehTLrCbJc/JDhNwC68k9336gVrAABvXQrRCqgZzNiO/UP9QlIND5BhlEPT5ULfFVe8FHJ5x0I0ikvrCcC+pyXMsJqK3trVHIj+xRZDJCAXv55I6CXd47EBnpV6j2U+/sNqFNu5RZCoURK/feQWmuyvAWSQXQE9NIlQUqgl25W5KmsKUA79LRDX+LWPAG9/FNGQC/vHIkJ9IqsdFQ4Phjq1DtWRUh4ZgkyqvSTt0jFHB0BvXQJIKCXTmvyJJ4CBPQE9OKtLpEsE9CLJKyAZgnoBRRTBFNiAr3bvVXwvPZ5gVEbtRXxuNdpEWZWekwS0EuXSwJ66bQmT+IpQEBPQC/e6hLJMgG9SMIKaJaAXkAxRTAlJtB7n38LLlG7bUb9uOdpGJ0r2uxXVjsQ0EuXeQJ66bQmT+IpQEBPQC/e6hLJMgG9SMIKaJaAXkAxRTDlmZ4M97s3kNygKVLVLoJ68D3cE+onoTZtxrdZicyKXW32K6sdCOilyzwBvXRakyfxFCCgJ6AXb3WJZJmAXiRhBTRLQC+gmCKY8nTVwN1FjeQ0PVLTswT14HHjS7iH/GjDpgKP+l6FSeCLCUEnUszGCOilSwABvXRakyfxFCCgJ6AXb3WJZJmAXiRhBTRLQC+gmCKYEhPoNYkXUeH48wVGnVGpBxJaLxNhZqXHJAG9dLkkoJdOa/IkngIE9AT04q0ukSwT0IskrIBmCegFFFMEU2ICPQvX6+J7cH2wzWLkJpUrorsdhNG5sggzKz0mCeilyyUBvXRakyfxFCCgJ6AXb3WJZJmAXiRhBTRLQC+gmCKYEhvoYcyEx92lcL/7A2DUZc9A59MSSU3nIcujjgizKl0mCeilyycBvXRakyfxFCCgJ6AXb3WJZJmAXiRhBTRLQC+gmCKYEh3ozTEbM6FOfwCFPhUGFz8Ytb4izKZ0miSgly6vUgJ91uUYpP7nJLLuJMCYkAFVVQ9oe9SA24zWUHg4STdp8lTqFCCgJ6AvcYuagF7+KSOgl3eOJAN6ecsg6+gI6KVLjxRAb0rRIfWL00hffQUw5Z+bsqIrPD7rBO1ztaWbOHkqVQoQ0BPQl7gFTUAv/5QR0Ms7RwT08s4Pi46AXrocSQH0WTfiEN/jlwIn5dSpKrw3DSrSxMMiHuGLxRvwKDoezloNRg7ugYZ1a+Dz79bi42ljUb92dez6/RTq1qqGuPgkLFy5HSqlkn8+bcJQeLi74qslG+FfxRdjnu+JjEwdZn+zGi0a1caIQd1zxabPMmBD8EHsOXgaOr0eNav74aO3RmPNlv24ciOU23V1ccaM14Zj254/+WdZWQYM7d8ZQ/p1wp17kZi78GcoFQpu9/n+nTGgZ3us2LgHDx7G8HjVKhUs+flwyigc/+sKNu04jPRMHfwqlccnb4+Hf+UKFvWLehSLj79cjncmDUeTBoFIS8/E5t8OY/hz3eDqos0eY6nfnPlr8TA6Dgoo0KtLKx57ekYm/v3NGkQ9jsWHU0YjM1OPn9bvxLyZk1DZtxxW/bIXz/XpiHLeHth3+C+s3XqA+/n6X68jIjIa8xZvQEpqGpo3qo2ZU8dwnYraCOgJ6Iu6hiQfT0AvueSFdkhAX2jJJB1AQC+p3A45I6B3SDaHBkkB9GmL/kbqPNtvR/a9MQEKr38AszATYoA4Z8E6TH11KAdbg8GI0PAo6HRZ2Lr7KLw83DDllSHYtusoGtcPRHRcAjffrUMLbAg+BE93Vwzo1R5fLtmI6NgEfPTWGA7WS9fuQOvmDfDqqGdzhbPn0BkOtOxzhUKB2PgkKJVK/LztAEYN7gHf8t7Z/Reu2MY/YxcMc+avw+RxA5GckoarN0NzXSgkJKVgyapfOfi/OvpZVPOrCEt+wiOjsWP/ccyaPh4atQrpGTqYTEarYLx9z59gYxhUTxo7EGnpGVi2fjcmjumfa0xB/VyctRzOb9y5j7cnDcP12/ez4z988gL+PH0ZgdWrYNwLvWGeb+SjWJz++zomjhkAlUrJ9WAXKOwihv2bXdAM6hOERvVqFCbVFvsS0BPQF3kRSW2AgF5qxQvvj4C+8JpJOYKAXkq1HfNFQO+Ybo6MkgLok6cdQsaWmzbD89k5FJpnHDsB6tqtMA65M6eOzeWHfX75RgiiYxMR1KYJ7t57kAvoO7dtlg3tbVrU5zBatYovMjL1SEhMhm8FH76bnBforcEouxuQmpYOF60WWq0G44f1wbbdRznQq1QqfPvjZrz72ggw2GUA3aZFAw7arZvXx9mLt/DwcSzc3VzxJD0dwwZ0tQi9iUmpmP7pYnTt0BxBrZsgoFplDvaWGrvL8N1PWzByUHe++8/uGGidNPmA3p5+T9Iy8J8F6/DO5OF4HJOQC+gzMnT46+JNvDyiL37dd5zPl11IsR35+w8e81gZ7DtrnZCU/ASrN+/jdzEmj3/Oauw2F0yODgT0BPSFWS+y6EtAL4s0FBgEAb28c0RAL+/8sOgI6KXLkRRAn/L+EaT/fM3mpMQCerYT3r5VI6zffpDverdoXIfv0G8MPsTLT7p1bIFXRvbjO+0M6Af37YQFy7fykpCqfr4Ivf8Qg/sG4Z3ZS3g5z+wZL4HtSlvaXZ6/bCv6dmuD8j6e3J6Xpzu+WLQeSSlPEBIWic/ee4WXvbALjTMXrmNgrw58Z5/dQfh+za9o17IR3N2c+V2DWdPH8QsAS34ydXpcuhaCk+eu4syFG5jzwauoFeCXT+M79x4geM8xvDyyHzb+eggtGtfFM03r5AN6e/qxnX12h4Ht0LMLJPMdBqYFa57ubjh/5TYY+LOSpZ9+3olne7TjOjLtnTRqXubDdulj4xKx8pc9GPFcd9Su6W9zbdjqQEBPQG9rjcjuewJ62aUkX0AE9PLOEQG9vPNDQC9tfqQA+vS1V5Hy4dGCJ6YAKt6dDLioHRKAlbzMW7Se79AzmGaN7QSzshkzeC5bvwvnr9zBlJeHZJfctGpaj9fJvzS8Lwdtc7kIqxP38nDH+au3OdDn3aFncMxKR8y19TqdHgajCcvW77RacsMuBtjONLsgYCUwOUtuWJyLVmxH96CWPPYjJy9i7NBe/O5CXj+paZlwd9XCyUnD+7J5secAOrVtmk+75Rt28wsLdteBaRQe+RhTX30eyzfsyVVyY0+/y9dDELzvOK/vv3U3Ih/QM////f4XPrfP3nsZe/84w58t6NS2CX7/8xwSElM40JtbQXEXdhEQ0BPQF3bNFHt/AvpiT4HNAAjobUpUrB0I6ItVfruc0w69XTIJ0kkKoDclZiKuywYYY9KsxuwyoSk/6aYo7diZy1i8Mpg/9BoTl8Br39u1bJgNng+j4/ku+6xp43LV0LOHQdmDmh9NGc3LRHLWwLPdZ0tAb67ZZzv8Hm4uiE9Mxsyp47Bm8z6rJTesrv7A0XM4d+km+vdoj5t372dfELBaeRNM6N+jHZfg2JkruH47DKOH9ODPBuT0M35YX8xftgXV/Svxvqy+/f03RvIa/ZyNxcjuGEydMJTfAWBlNaxUiJW+LF37Gx/HHryt5OuDx7EJfOfdUj8njQb3wqMQGOCHD94cDYVSga9/+AX3wh/xh1ofxcRzt+x5BHbngV0gfT9vOvT6LMxduB6eHq5IeZLOLwSOnryIE+euQq1iF24m/Pvdl/PF7cgaKNFAz259TJm5gN9qMbfFc6dxQc2N1Uux2zWstW3RAOx789PEUXHpFjVjP9zWvnNEZBojrAIE9MLqKYY1AnoxVBXOJgG9cFqKZYmAXixl89uVAuiZV92fEUgc9ZvFYyvVzSqi3G9DAc3TByeL2hj8ujg78Z1psRsrH4HJBI3GsTsL9sZnyQ8ru2GN1cSz75OSU2Ey/XMuKCvl8fZ0z34g1V5fQvZj8bB8sOcEzPlgsWZkZAoC8uZYSzTQx8Ql8gcdWI0Vg3R2FTlz7jIs/+Z9/sQw+/eXizdi3aKZ/GnrmfOW8XnP/Wgi/5OAXsglK50tAnrptHbUEwG9o8pJM46AXhqdi+KFgL4o6hVurFRAz6Ji59Gn/XABmYfDYQhLglPbKnAe3Qja3kU/5aRwsy59veUK9FIpXaKBPq9IDPDHvTUXH0wZxXfpGcDXquGfXfeVF/AJ6KVaZsL6IaAXVk8xrBHQi6GqcDYJ6IXTUixLBPRiKZvfrpRAL92syFNZU6BUAT2rW2K1Yd/OfhM1q1fm5Tgd2zTJBvqc37Md/JjEDIv59vV2tvpdWVsgcpyvk0YFZyclkp88vdVWJpr4d00FlZEBY4bOCB27FSvjxm6FSnFLWm4SuDmr4apV40lGFtIys+QWHsUDwNlJzasvUjL+lx8LbxcloYRRgP2fT40UKOkKlBqgN9fTmwHe/O9xw3pn19TnBXpdltFi/pzUSlj7rqQnvDTEr1SAP/GeZSxD/8OVsKmqVQoYTSYYLf+IyWYZZhmMUP/vZR+yCUqCQFQqBVRKBQwGEz+Vgpr8FFAqwS82WY54K2EX9fJT1HpE7P98aqRASVdAtkDPymV27D9hUV92oH/O45PM8F65Yrns+vi8gM8M5QV6KrkpmcuXSm7knzcquZF3jqjkRt75YdFRyY10OaKSG+m0Jk/iKSBboLd3ypZg3jyWaujtVbFk9SOgl3++COjlnSMCennnh4Be2vwQ0EurN3kTR4ESDfSWduFzykSn3IizaIrbKgF9cWfAtn8CetsaFWcPAvriVN8+37RDb59OQvQioBdCRbJR3AqUaKBnJTQTZnyF5NTcL2oY1KdjdukNnUNf3EtMeP8E9MJrKrRFAnqhFRXWHgG9sHqKYY2AXgxVLdskoJdOa/IkngIlGujFk4UskwKkAClACpACpAApQAqQAiVDAQL6kpEnipIUIAVIAVKAFCAFSAFSgBSwqAABPS0MUoAUIAVIAVKAFCAFJFAgFPFYZ7iIB6YkpCATvnDDM0p/DFc2gSs0EkRALkqrAgT0pTWzNC9SgBQgBUgBUoAUkIUCadBjg/ES9hlvw9KbH7zhjFdUrdBRUV0W8VIQJU8BAnobOct7Hv7iudOyX1RV8tJd+iJmJxmxNwKbW9sWDcBy5OpCb/6TW7bZA+ohYZHZD6zLLb6yEk9MXCLGvTUXEVHRfMr0O02emWeHPsz+ZjW+nzcdvuW95RkkRWW3AveRiHey9hTYv6miMj5VdbfbZt6O98IfYs6CdfzjWdPG4VF0PBau3A4njRotGtfBpLEDcOHqXfywZgc/TMTDzQUzp46FRqPG7ZAIDOjVvkDfJ85etTh24fJtmPX2OFT2LYfT56/jp3U78d9/vc5tsta+VSP+56pf9uK5Ph2xdddRHDh6Fr27tMbkcQNz+fx66SZcuRGKrCwDhvbvjCH9OvEXrIVFPMIXizfwOTlrNRg5uAcG9w1C1KNYfPzlcrwzaTiaNAjMtpWSmoYvl2zE+GF9UDewqsOalqSBBPQFZIv9x/fdT1swa/o4Doh5X0xVkhJdWmNlkBgYUIVfZBX0ToLSOv+SMK+cF105T6AqCbGXthjzHvVLv9Pkl+GcF1zV/Cpi3aKZBPTyS1OhI9puvIb1xks2x61VvwA3ONnsZ60D+33LGvs/0fz3ru2bY0PwIXi4u+DkuWv44M1R8PHyQEamDgaDkcPy1ZuhGDHI+sVEQlIKB2RLYxevCkb3oBYcvr9Zuhm3QyPwxcxJuHrrXnYs7C8LV2zDqME94ObqjGXrd2PimP75Nt/MfTzcXTFn/joO/N6e7vxCZeqrQ+FfuQKPOTQ8CnVqVsX2PX8iPDIari5aTBr79OJAn2XAklXBeByTgLFDe6FRvRoO61mSBhLQFyJb5l+0H0wZRbv0hdBNyq4M8E/8dYV26aUU3U5ftENvp1Aidsu762vrXR4ihkKmbShAO/Sla4ksMpzCEdNTwC2ozVP1Rl1FBVvdrH5vDeiXrPoVdQKr4uCxc/Cv7IseQS1Rt1Y1aJ00fLPSFtBn6vSY9eVyi2PPXrqJh4/j0KpZPQ7R12+H4e1JwzjQbww+BL9KT+cT8TCag749QK9SqfDtj5vx7msjEPkoFjv2H+d3E3I2dkHCNl1HDuoO9v/LjNeGc/jfuvsoGtQJwLWb99C4fiABvcOrqRQPZD8oM+cuw/Jv3i8zC6SkpZOVSLE296OJJS30Uh8vAX3xpzjvy/ZYRPQzU/x5sRQBAb088+JoVD8a/8IB412bw4UG+gXLtsHdzQU9Oz+DMc/34v5vhYTjrws3cODIWUyb8AI8PdxsAj0bx3a+rY31r+KLlRv34ONpYzlcm4E+OSUNHXKU3Lw8sl8uoL90LYSXllWuWA7fzn4T36/+FUkpT3h55mfvvcLLaNjPgiWgv3PvAYL3HAOzufHXQ2jRuC7fqWe79+yuxI79J9Cobg20adEAKpXSpvYlvQPt0NuRwZwvsKJ6UzsEK6YulmClmEIhtxYUIKAv/mXBfkbWbTmQ6w4WAX3x54WAXp45EDKq/cY7+Ml4tkCTCgDr1SOghcph15Z26Fn5jRnI9Xp9dpnL73+eQxTbWW9azybQM5gvaOygvkG4/+Axqvn58lIZM9Azv2b/hSm5YbXyqzfvw+wZL4HdHZi3aD3foS/v48nnkpT8BFt2HeH19VWr+CI2PgnhkY/x8oh+uHQ9hPc5ee4q6gZW47X4zlrHy5gcTobEA8sk0JtvM5+5cMOi3NagnUpupFud1t4CzCKwVFdKd0+ky43ZU94HxnNG8M7k4Xh11LO5giKglz5HeT3SDn3x58DeCGiH3l6l/q+9Ow+xqgzjOP5oUTqa7RJZYSlF2YLtZCvSQhFBZWEkrSa0mCllTRG2TRhUVhqVbaBlZVBitAjZQoZFZRlm0UIRWRljJjZZVMZz4r289533nPfMnTtz7zvv9/4j3LnnnPf9POfK77zznDNxfG6j/CVX//2SrJdNuQM+rf8+cnH/Q2qekLamzJzzdLb99CvOy3rZ7UC9/reNMu2WB6WlZUB2Q6zeOHrjlAnSvm6DzJr7fNafrq+xxxwsxxxxYNU4ira123U0X4UC/dJ3PpLnFr8p55x+fKe+fRP69UbwJW99IB988rlce/l4Wf7hKpn9+AtZm9Av7b/KAfuOkHW/bpDJl54l224zKLsfoO3++VkfvQZ8fT27aCktNzWfTQlsqCFmxPBhncJKAlNv2ikS5pu2NAT6JisNPfRNVpCC4RDo46lV2ZGu3PyT3PrPUu9jK0f220Hu2OIk2VJ6vjXk/9X2v7P2FN9LQ/nvHdUXHltvvZUMGdyStd0UbVvWotbPdfzxpwwcsFW2Ms+rWiDJFfqyJ4Hbt2VWjdtaJ3JTbFnEHv4cbTY9DFzH3bNCX0fMGnfFU25qhGvAZgT6BqD3wiH1efSL/l0tH29eIz9u3ij79dtZxvYfKYf1G9YLRy93iKJAX24PfKoRAgT6AnVfaw499I04TfOP6Wv70FUEblxunjq5fytAR8b3qHH14Tn0jbMvc2S3ProNj3stI8dnEEhbgECfdv2ZPQIIIIAAAggggEDkAgT6yAvI8BFAAAEEEEAAAQTSFiDQp11/Zo8AAggggAACCCAQuQCBPvICMnwEEEAAAQQQQACBtAUI9GnXn9kjgAACCCCAAAIIRC5AoI+8gAwfAQQQQAABBBBAIG0BAn3a9Wf2CCCAAAIIIIAAApELEOgjLyDDRwABBBBAAAEEEEhbgECfdv2ZPQIIIIAAAggggEDkAgT6yAvI8BFAAAEEEEAAAQTSFiDQp11/Zo8AAggggAACCCAQuQCBPvICMnwEEEAAAQQQQACBtAUI9GnXn9kjgAACCCCAAAIIRC5AoI+8gAwfAQQQQAABBBBAIG0BAn3a9Wf2CCCAAAIIIIAAApELEOgjLyDDRwABBBBAAAEEEEhbgECfdv2ZPQIINJnAYwtelmXvfyqz266WloEDGjq6N95dIVe23idTJ50jl4w/tS5j6fhjU7bPNT+3y7wHWmXnHbery37ZCQIIIJCyAIE+5eoz94YJrPriW7l02l2y7ZDBnUJN0c96e8AmfI05/IC6BbrenkPZ42l4nTl7QcNDZijQm5rsMnQHabthYtnpdflzv7SvlwlXtcm404+vqr15//s1a7N96oXHCUeNrtq/uRCw39x916EV296aQ5cnzQYIIIBApAIE+kgLx7DjFjChfcPGjk6BqPXOubLotWViB6BGzZZA3/vyoUDfWyPyjcM9HzTcX37DLJkx7UIZtc/wbGh6/r6xbIU8evd1lff0fd3f19/+ULkI0e/A1Blz5J4ZV1R9rrfmx3EQQACBviRAoO9L1WQu0QiYMHPq2CPlk1VfVdor9P0Zdz8po/cfKW8vX1m1WmxfBOhEzzh5TNUKrbkQMAi+FVFdaddQpRcM+jpi9L6FrR3uPocMbqkENXel1t6X+dlF554iTzz7qpjVXG3dOPLg/bLfTujFjLvCa0LkmacdK9Nvfzj7ue/CRj93z8PPVeptrxKblfbpV46X1ra52XF8x9WNTSuJa2t8b5oyIWsPcX9DoS760hVyO+QaW3vMRWP1nbDG4KBRI+WR+Ys7Gfgusopqrzsw27y3YnXlkEVtNHkXcub8fPDOKZVWGT32iccdmq3Sq72au2E+74tpO0bz5WWgCCCAQBMKEOibsCgMqe8LmEB//VXnZW0eGj41EJlVzBHDh8nCxW9WAr27mulrWWi7f76ccfLRVSulP61dlwV2fWkw1UBnwm9eS4WtnxfszLZm3LqNhjNzvN87NmXtGvoyfdKmDcMOu26biwm/9sWKvV/tKXdXjk0Yb2udWAmVOlf3YkU/t+i1d6R18vnZuHzbuS03efP3BXrb1hiGxpoX6PVixQ7ctoGppX2RUVR7NXODs9bv3kcWil6w+Pr081bPQyv0bq1C32St/7yFS5rifoHQWPk5Aggg0MwCBPpmrg5j67MCdmBa/tFn2U2QN0+9QK697aGsfUHfswO9BiUN+faNiaGeb3s1dVDLgOBKsw87L9C67RMmIJsWiqE7bZcFejvw+y4C3Pd8bR62le7XbfEwFxP6r66Yh1zsedquvu26EujdVXxfO4o71rxA794Uaxvsuccu3lra+ypT+6IvV1HQzuuhr6Uv3rfi32e/9EwMAQQQ6EEBAn0P4rJrBPIE3JCq4XfwoIGy9167ZaFUg60J9CaM2+0SZr++thr7c6ZFJi8EhloeigKtadux52iOV89Ab4f+oTtuX9WuYx/brOoXBXrfzZpF23Un0PvaeMx43XYpex6+ixrbQH/z4LYB+Vpq7PYou+0n1GalY6ll5byWQJ930cP/HAgggAACXRMg0HfNi08jUBcBt6VBA9ejT71U6T32BfqiJ82Y8HjCmNHemw57ItCbFXEfSJnVeN2uzAq9G+hDN1LmBXrfzZr2BU29V+hrvemzq4E+VHtzs6p6lw32tQR63X9XW24I9HX574SdIIAAAkKg5yRAoAECbtjTYLPgxddl8iVnZaOxA70+pzu0ku4LYGXaNEL7NSHNbfcJPYmlnoHe99sMu5XHLV9RMJ8w7qSqRyza8/e1f3Rnhd5nUOZUC7UduRdnodrbgd4cP9SWVGugL7op1j3HdSy03JQ5I/gMAgggEBYg0IeN+AQCdRcIrd66gd73B37sGxu1zcZ+uogJoqu//C5b9a91hT5v1dW3KqzHvH3WPLnmsnGZVz166H1tHL6VdvX55rsfs3sMigK9/ex2Y2raX9ybZE3R3VVnd7uiR3uGxuo7sdxA7xq4x3NDtFt7bX/SJwbNvGlS5ck0oQuy0PlZ9IXwzdl3vug+ar1wqPsXkh0igAACkQsQ6CMvIMOPUyAUmNxAb1Yz7cc96nv2k1Dsdgrtn54y8Wx5/JlXsud8dyfQ2zdBFj22UsdjwnF3Vujtx1Ha+7Qr7T4K0h5X3uqzezOnjtW8zB9osvdr5uL2p7vbhZ7VXzTWvEBfZOA7XlHtdYU+9FhLdxyhOYW+db57FXyPyfTd7B3aNz9HAAEEEOgsQKDnrEAAgaYRCK0cN81AExhIT9cidFGbADFTRAABBOomQKCvGyU7QgCB7gr0dIjs7vhS2r7M3ymo1aOWJ+LUeiy2QwABBFIQINCnUGXmiEAkAgT65iqU796N7o7QhPk1P7dX/SXk7u6X7RFAAIGUBQj0KVefuSOAAAIIIIAAAghEL0Cgj76ETAABBBBAAAEEEEAgZQECfcrVZ+4IIIAAAggggAAC0QsQ6KMvIRNAAAEEEEAAAQQQSFmAQJ9y9Zk7AggggAACCCCAQPQCBProS8gEEEAAAQQQQAABBFIWINCnXH3mjgACCCCAAAIIIBC9AIE++hIyAQQQQAABBBBAAIGUBQj0KVefuSOAAAIIIIAAAghEL/Af8FDctLaIWSAAAAAASUVORK5CYII=",
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
    "server = app.server\n",
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
    "    app.run_server(debug=True)"
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
