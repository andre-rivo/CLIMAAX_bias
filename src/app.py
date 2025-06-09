#!/usr/bin/env python
# coding: utf-8

# In[ ]:





# In[1]:


import numpy as np
import geopandas as gpd
import os
import xarray as xr
import matplotlib.pyplot as plt
import glob
import rioxarray as rio
from tqdm import tqdm
import plotly.graph_objects as go
import pandas as pd
import plotly.express as px
import plotly.colors as pcolors
from dash import dcc, html, Input, Output, Dash


# In[3]:


#import eobs tas

coords = ['longitude','latitude']
eobs_tas = xr.open_dataset('tg_ens_mean_0.1deg_reg_v30.0e_REMAP_FAIZAN_mean1971-2000.nc')
eobs_tas=eobs_tas.assign_coords({coord: eobs_tas[coord].round(3) for coord in coords})
    


# In[4]:


#import eobs pr
coords = ['longitude','latitude']
eobs_pr = xr.open_dataset('rr_ens_mean_0.1deg_reg_v30.0e_REMAP_FAIZAN_MEAN_YEARSUM_1971-2000.nc')
eobs_pr=eobs_pr.assign_coords({coord: eobs_pr[coord].round(3) for coord in coords})
    


# In[5]:


gcm_list = ['EC-EARTH','HadGEM2','NCC','CNRM-CERFACS','IPSL','MPI']
rcm_list = ['CLMcom','SMHI-RCA4','KNMI','DMI-HIRHAM5','REMO2015','ALADIN63','REMO2009','WRF381P']


# In[6]:


def mean_bias_tas(string,eobs=eobs_tas,coords=['longitude','latitude']):
    model=xr.open_mfdataset(string)
    model = model.convert_calendar('noleap',align_on='year')
    try:
        model = model.reset_coords('height')
    except ValueError:
        print('no height')
    model = model.assign_coords({coord: model[coord].round(3) for coord in coords})
    mean_bias = model.tas[0] -273.15 - eobs.tg[0]
    return mean_bias


# In[7]:


def mean_bias_pr(string,eobs=eobs_pr,coords=['longitude','latitude']):
    model=xr.open_mfdataset(string)
    model = model.convert_calendar('noleap',align_on='year')
    model = model.assign_coords({coord: model[coord].round(3) for coord in coords})
    mean_bias = (model.pr[0]*86400 - eobs.rr[0]+0.00001)/(eobs.rr[0]+0.00001) *100
    return mean_bias


# In[8]:


#calculate tas bias in every grid cell
bias_list_tas=[]
#coords = ['rlon','rlat']

names_list = []

for gcm in gcm_list:
    for rcm in tqdm(rcm_list):
        try:
            string=f'models_faizan/all_projections/tas*{gcm}*{rcm}*1971_2000_timmean.nc4'
            files = glob.glob(string)
            bias_list_tas.append(mean_bias_tas(files))
            names_list.append(str(gcm)+'_'+str(rcm))
        except OSError:
            print(f"No files found for {gcm} and {rcm}, moving to next combination.")
            continue


# In[9]:


#calculate pr bias in every grid cell
bias_list_pr=[]
#coords = ['rlon','rlat']

for gcm in gcm_list:
    for rcm in tqdm(rcm_list):
        try:
            string=f'models_faizan/all_projections/pr*{gcm}*{rcm}*1971_2000_timmean.nc4'
            bias_list_pr.append(mean_bias_pr(string))    
        except OSError:
            print(f"No files found for {gcm} and {rcm}, moving to next combination.")
            continue


# In[10]:


#creates lists of bias for every NUTS2 region, skips region outside of CORDEX domain

regions_gdf = gpd.read_file('NUTS_2_SHAPE.shp', encoding='utf-8')
crs = regions_gdf.crs


#crs = ccrs.RotatedPole(pole_latitude=39.25, pole_longitude=-162)

#regions_gdf=regions_gdf.to_crs(crs)

regions_ex=regions_gdf.explode(column='NUTS_ID')

tas_bias_masked_list=[]
for model in np.arange(len(bias_list_tas)):
    for i in np.arange(len(regions_gdf)):
        try:
            bias_list_tas[model] = bias_list_tas[model].rio.write_crs(crs)
            cropped = bias_list_tas[model].rio.clip(geometries=[regions_ex['geometry'][i]], all_touched=True)
            tas_bias_masked_list.append(cropped)
        except Exception as e:
            print(f"Skipping region {i} due to error: {e}")
            continue
        #bias_median[i,model] = np.median(bias_masked)
        #bias_iqr[i,model] = np.quantile(bias_masked,0.75,method='midpoint') - np.quantile(bias_masked,0.25,method='midpoint')
        
pr_bias_masked_list=[]
for model in np.arange(len(bias_list_pr)):
    for i in np.arange(len(regions_gdf)):
        try:
            bias_list_pr[model] = bias_list_pr[model].rio.write_crs(crs)
            cropped = bias_list_pr[model].rio.clip(geometries=[regions_ex['geometry'][i]], all_touched=True)
            pr_bias_masked_list.append(cropped)
        except Exception as e:
            print(f"Skipping region {i} due to error: {e}")
            continue


# In[11]:


#calculate median bias inside each region

weights = [ np.cos(np.deg2rad(i.latitude)) for i in tas_bias_masked_list ]

combinations=20 ######## UPDATE this to become automatic i.e., calculate automatically how many gcm-rcm combinations exist


median_bias_tas = [i.weighted(weights[num]).quantile(0.5, dim=('latitude', 'longitude')).values for i,num in zip(tas_bias_masked_list,range(len(tas_bias_masked_list)))]
#median_bias_tas = np.array(median_bias_tas)
#median_bias_tas = median_bias_tas.reshape(combinations,-1)

median_bias_pr = [i.weighted(weights[num]).quantile(0.5, dim=('latitude', 'longitude')).values for i,num in zip(pr_bias_masked_list,range(len(pr_bias_masked_list)))]
#median_bias_pr = np.array(median_bias_pr)
#median_bias_pr = median_bias_pr.reshape(combinations,-1)


# In[ ]:





# ## Create map and attach plots to regions

# In[12]:


regions_gdf = gpd.read_file('NUTS_2_SHAPE.shp', encoding='utf-8')

regions_gdf = regions_gdf['geometry']

df=gpd.GeoSeries.get_coordinates(regions_gdf)


# In[13]:


#drop extra-continental regions
#regions_ex2= regions_ex.drop(index=[141,142,143,144,145,163,256,257]).reset_index(drop=True)

regions_ex2= regions_ex.drop(index=[141,142,143,144,145]).reset_index(drop=True)

#empty list to be filled with all plots produced
plots_list = []

#label for each model combination
label_list = names_list

cmap=plt.get_cmap('tab20')
color_list = cmap(np.linspace(0, 1,combinations+1))


#loop through all regions
for region in tqdm(range(len(regions_ex2))):


    #make one figure for each region
    fig,ax=plt.subplots(1,1,figsize=(6,4))
    #region_bias_tas = median_bias_tas[:median_bias_tas/combinations]
    #region_bias_pr = median_bias_pr[:,region]
    
    #in the same figure for each region, plot temp bias vs prec bias, with each color representing one model and having a label
    for model in range(combinations):
        index = region + ( int( len(median_bias_tas)/combinations) * model)
        region_plot=ax.plot(median_bias_tas[index],median_bias_pr[index],'o',color=color_list[model],label=label_list[model])    
        
    
    #Add 0,0 marker for observations
    ax.plot(0,0,'s',color='k',markersize=7,label='Observations')
    
    #Axes titles
    ax.set_ylabel('Mean precipitation bias (%)')
    ax.set_xlabel('Mean temperature bias (°C)')
    
    #Set origin to 0,0 and middle of the plot
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.plot(xlim, [0, 0], color='b', linewidth=1.5,label='0 Prec bias')
    ax.plot([0, 0], ylim, color='r', linewidth=1.5,label='0 Temp bias')
    
    
    #Add legend and append the plot for each region to the list
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.05, 0.95), loc='upper left', fontsize=6)
    
    
    plt.title(regions_ex2['NUTS_NAME'][region])
    
    fig.subplots_adjust(right=0.95)
    plt.tight_layout()
    
    plots_list.append(fig)
    plt.close(fig)


# In[15]:


#### PLOTLY

label_list = names_list

# Drop extra-continental regions
regions_ex2 = regions_ex.drop(index=[141, 142, 143, 144, 145]).reset_index(drop=True)

# Empty list to store plotly figures
plots_list_plotly = []

# Color list for each model combination
colorscale = pcolors.qualitative.Light24
color_list = [pcolors.sample_colorscale(colorscale, i / combinations) for i in range(combinations + 1)]

# Loop through all regions
for region in tqdm(range(len(regions_ex2))):
    region_name = regions_ex2['NUTS_NAME'][region]
    
    # Data collection for plotting
    x_data = []
    y_data = []
    model_labels = []
    colors = []

    for model in range(combinations):
        index = region + (int(len(median_bias_tas) / combinations) * model)
        x_data.append(median_bias_tas[index])
        y_data.append(median_bias_pr[index])
        model_labels.append(label_list[model])
        colors.append(color_list[model])

    # Create Plotly figure
    fig = go.Figure()

    # Add model data points
    for i in range(combinations):
        fig.add_trace(go.Scatter(
            x=[x_data[i]],
            y=[y_data[i]],
            mode='markers',
            marker=dict(color=colors[i], size=10),
            name=model_labels[i],
            customdata=[[region, model_labels[i]]],  # For interaction
            hovertemplate=(
                f"<b>Model:</b> {model_labels[i]}<br>"
                f"<b>Temp Bias:</b> {str(x_data[i])[:4]} °C<br>"
                f"<b>Prec Bias:</b> {str(y_data[i])[:4]} %<extra></extra>"
            )
        ))

    # Add origin point (observations)
    fig.add_trace(go.Scatter(
        x=[0], y=[0],
        mode='markers',
        marker=dict(color='black', size=10, symbol='square'),
        name='Observations',
        hoverinfo='skip'
    ))

    # Add reference lines at zero
    fig.add_trace(go.Scatter(x=[min(x_data + [0]), max(x_data + [0])], y=[0, 0],
    mode='lines',
    line=dict(color='blue', width=1.5),
    name='0 Prec bias',
    showlegend=True))


    fig.add_trace(go.Scatter(x=[0, 0], y=[min(y_data + [0]), max(y_data + [0])],
    mode='lines',
    line=dict(color='red', width=1.5),
    name='0 Temp bias',
    showlegend=True))

    # Layout
    fig.update_layout(
        title=f"{region_name}",
        xaxis_title="Mean temperature bias (°C)",
        yaxis_title="Mean precipitation bias (%)",
        legend=dict(font=dict(size=8), x=1.02, y=0.95),
        margin=dict(l=40, r=40, t=40, b=40),
        height=400,
        width=600)

    
    # Append to list
    plots_list_plotly.append(fig)   
    


# In[66]:


plots_list_plotly[100]


# In[16]:


# Initialize Dash app
app = Dash(__name__)
server = app.server

# App layout (select region and show main graph + details)
app.layout = html.Div([
    html.H4("EURO-CORDEX Climate Model Biases and Uncertainty at NUTS2 level"),
    dcc.Dropdown(
        id='region-select',
        options=[{'label': name, 'value': i} for i, name in enumerate(regions_ex2['NUTS_NAME'])],
        value=0
    ),
    dcc.Graph(id='main-graph'),
    html.Div(id='detail-output')
])

# Main graph update on region selection
@app.callback(
    Output('main-graph', 'figure'),
    Input('region-select', 'value')
)
def update_main_plot(region):
    region_name = regions_ex2['NUTS_NAME'][region]
    x_data, y_data, model_labels, colors = [], [], [], []

    for model in range(combinations):
        index = region + (len(regions_ex2) * model)
        x_data.append(median_bias_tas[index])
        y_data.append(median_bias_pr[index])
        model_labels.append(label_list[model])
        colors.append(color_list[model])

    fig = go.Figure()

    for i in range(combinations):
        fig.add_trace(go.Scatter(
            x=[x_data[i]],
            y=[y_data[i]],
            mode='markers',
            marker=dict(color=colors[i], size=10),
            name=model_labels[i],
            customdata=[[region, model_labels[i]]],
            hovertemplate=(
                f"<b>Model:</b> {model_labels[i]}<br>"
                f"<b>Temp Bias:</b> {str(x_data[i])[:4]} °C<br>"
                f"<b>Prec Bias:</b> {str(y_data[i])[:4]} %<extra></extra>"
            )
        ))

    fig.add_trace(go.Scatter(
        x=[0], y=[0],
        mode='markers',
        marker=dict(color='black', size=10, symbol='square'),
        name='Observations',
        hoverinfo='skip'
    ))

    fig.add_trace(go.Scatter(x=[min(x_data + [0]), max(x_data + [0])], y=[0, 0],
                             mode='lines',
                             line=dict(color='blue', width=1.5),
                             name='0 Prec bias'))

    fig.add_trace(go.Scatter(x=[0, 0], y=[min(y_data + [0]), max(y_data + [0])],
                             mode='lines',
                             line=dict(color='red', width=1.5),
                             name='0 Temp bias'))

    fig.update_layout(
        title=f"{region_name}",
        xaxis_title="Mean temperature bias (°C)",
        yaxis_title="Mean precipitation bias (%)",
        legend=dict(font=dict(size=8)),
        height=400,
        width=600
    )

    return fig

# Uncertainty plot when clicking a model point
@app.callback(
    Output('detail-output', 'children'),
    Input('main-graph', 'clickData')
)
def display_details(clickData):
    if clickData is None:
        return html.Div("Click a model to see more details.")

    model = clickData['points'][0]['customdata'][1]

    detailed_fig = go.Figure()
    detailed_fig.add_trace(go.Scatter(
        x=np.arange(10),
        y=np.random.rand(10),
        mode='lines+markers',
        name=f'{model} Time Series'
    ))
    detailed_fig.update_layout(title=f"Uncertainty for {model}")

    return dcc.Graph(figure=detailed_fig)

# Run app
if __name__ == '__main__':
    app.run_server(debug=True)


# In[ ]:


























