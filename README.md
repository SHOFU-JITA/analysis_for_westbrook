# analysis_for_westbrook

This repository contains the analysis scripts.
---

## 1. Overview & Reproducibility

This workflow allows users to reproduce all figures, tables, and statistical results presented in the paper. 

- **Primary language:** R (version 4.6.1) / Python (version 3.10)

## 2. Directory Structure

```text
├── CJS/
│   ├── CJS_bkt3.R   #CJS analysis, Generates Figure 1b,c
│   ├── CJS_bkt4.R   #CJS analysis, Generates Figure S3b,c
│   └── CJS_bkt9.R   #CJS analysis, Generates Figure S4b,c       
├── ClusterChange/
│   ├── Clusterchange_bkt3_rev.R  # Analysis for Temporal Shifts in Growth Clusters, Generates Figure S6 
│   └── Clusterchange_bkt4and9.R  # Analysis for Temporal Shifts in Growth Clusters, Generates Figure S5
├── ClusterChange_Env/
│   └── logit_temp_flow_density.R  # Analysis for Effects of Density, Flow, and Temperature on Growth Clusters, Generates Figure 3c, Figure S7, Table S7
├── Clustering/
│   └── Clustering_bkt.R  # Growth cluster analysis, Generates Figure 1a, Figures S3a,S4a, Tables S5,S6
├── Combining/
│   └── Combining_bkt_rev2.R  # Raw data filtering and formatting, combining clustered data and PIT-tag antenna data
├── CV_window/
│   └── CV_combined_bkt3_tail.R  # Analysis for Life-history variation and portfolio effects through time, Generates Figure 4
├── Environment/
│   └── Temperature 
│       └── Temperature_seasonal_rev.R  # Analysis of water temperature, Generates Figure S1 and Tables S1,S2
├── Flow/
│   ├── analysis
│   │   └── Flow_analysis.R  # Analysis of flow, Generates Figure S2 and Tables S3
│   └── hindcast
│       └── hindcast_rev.R  # Analysis for hindcasting flow
├── JS/
│   └── JS_bkt.R   # Analysis for Population Dynamics, Generates Figure S8  
├── Movement/
│   ├── Distribution_bkt3.R   # Raw data filtering and formatting, Generates Table S4
│   └── Distribution_bkt4.R   # Raw data filtering and formatting, Generates Table S4
├── Movingout/
│   ├── Movingout_bkt3.R   # Raw data filtering and formatting for CJS
│   ├── Movingout_bkt4.R   # Raw data filtering and formatting for CJS
│   └── Movingout_bkt9.R   # Raw data filtering and formatting for CJS
├── PopulationDynamics_cluster/
│   └── PopulationDynamics_cluster_graph_bkt3.R  # Analysis for Temporal Shifts in Growth Clusters, Generates Figures 3a,b  
└── README.md
