import pandas as pd
import json
import numpy as np

def generate_html_report(df, output_path="report.html"):
    """Generates a comprehensive HTML report from the results DataFrame with interactive charts and DataTables."""
    
    # --- PREPARE DATA FOR CHARTS ---
    
    # 1. bRo5 Status Distribution
    if 'bRo5_Status' in df.columns:
        status_counts = df['bRo5_Status'].value_counts().to_dict()
    else:
        status_counts = {}

    # 1b. Lipinski Status Distribution
    if 'Lipinski_Status' in df.columns:
        status_counts_lipinski = df['Lipinski_Status'].value_counts().to_dict()
    else:
        status_counts_lipinski = {}
        
    # 2. Chemical Space (MW vs LogP)
    scatter_data = []
    if 'MW' in df.columns and 'LogP' in df.columns:
        for _, row in df.iterrows():
            # Use 'Compound' name if available, else SMILES
            label = str(row['Compound']) if 'Compound' in row and pd.notna(row['Compound']) else row.get('SMILES', '')[:20] + '...'
            
            scatter_data.append({
                'x': row['MW'],
                'y': row['LogP'],
                'name': label,
                'status': row.get('bRo5_Status', 'Unknown'),
                'lipinski': row.get('Lipinski_Status', 'Unknown')
            })

    # 3. Similarity Analysis (Hybrid vs Tanimoto vs LLM)
    # Find columns related to scores
    hybrid_cols = [c for c in df.columns if c.startswith('Hybrid_')]
    tanimoto_cols = [c for c in df.columns if c.startswith('Tanimoto_')]
    llm_cols = [c for c in df.columns if c.startswith('LLM_Cos_')]
    
    sim_data = []
    if hybrid_cols:
        # Use the max score across all controls for each compound for visualization
        df['Max_Hybrid'] = df[hybrid_cols].max(axis=1)
        df['Max_Tanimoto'] = df[tanimoto_cols].max(axis=1) if tanimoto_cols else 0
        df['Max_LLM'] = df[llm_cols].max(axis=1) if llm_cols else 0
        
        for _, row in df.iterrows():
            sim_data.append({
                'hybrid': row['Max_Hybrid'],
                'tanimoto': row['Max_Tanimoto'],
                'llm': row['Max_LLM']
            })

    # 4. Docking Data
    docking_data = []
    docking_table_data = []
    top_docking_candidates = []
    has_docking = 'Docking_Score' in df.columns and df['Docking_Score'].notna().any()
    has_vina = 'Docking_Score_vina' in df.columns
    has_autodock = 'Docking_Score_autodock' in df.columns
    
    if has_docking:
        # Get docking stats
        docking_scores = df['Docking_Score'].dropna().tolist()
        min_dock = min(docking_scores) if docking_scores else 0
        mean_dock = sum(docking_scores) / len(docking_scores) if docking_scores else 0
        
        # Prepare detailed docking table data with compound names
        docked_df = df[df['Docking_Score'].notna()].copy()
        docked_df = docked_df.sort_values('Docking_Score')
        
        for idx, row in docked_df.iterrows():
            # Use 'Compound' name if available, else SMILES
            compound_name = str(row['Compound']) if 'Compound' in row and pd.notna(row['Compound']) else 'Unknown'
            smiles = row.get('SMILES', '')[:40] + ('...' if len(row.get('SMILES', '')) > 40 else '')
            
            entry = {
                'name': compound_name,
                'smiles': smiles,
                'hybrid': round(row.get('Max_Hybrid', 0), 3),
                'docking': round(row['Docking_Score'], 2)
            }
            
            # Add mode-specific scores if available
            if has_vina:
                entry['vina'] = round(row['Docking_Score_vina'], 2) if pd.notna(row.get('Docking_Score_vina')) else 'N/A'
            if has_autodock:
                entry['autodock'] = round(row['Docking_Score_autodock'], 2) if pd.notna(row.get('Docking_Score_autodock')) else 'N/A'
            
            docking_table_data.append(entry)
        
        # Top 10 docking candidates
        top_docking_candidates = docking_table_data[:10]
        
        # Prepare scatter data (Hybrid vs Docking)
        for _, row in docked_df.iterrows():
            label = str(row['Compound']) if 'Compound' in row and pd.notna(row['Compound']) else row.get('SMILES', '')[:20] + '...'
            
            docking_data.append({
                'hybrid': row.get('Max_Hybrid', 0),
                'docking': row['Docking_Score'],
                'name': label
            })

    # --- HTML STRUCTURE ---
    html_string = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Comprehensive Analysis Report</title>
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width, initial-scale=1">
        
        <!-- Plotly.js -->
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        
        <!-- jQuery & DataTables -->
        <script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>
        <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
        <script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>
        
        <style>
            body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; padding: 0; background-color: #f8f9fa; color: #333; }}
            .header {{ background: linear-gradient(135deg, #2c3e50, #3498db); color: white; padding: 20px 40px; margin-bottom: 30px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }}
            .header h1 {{ margin: 0; font-size: 24px; }}
            .header p {{ margin: 5px 0 0; opacity: 0.8; }}
            
            .container {{ max-width: 1400px; margin: 0 auto; padding: 0 20px; }}
            
            .card {{ background: white; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.05); padding: 20px; margin-bottom: 20px; }}
            .card h2 {{ margin-top: 0; font-size: 18px; color: #2c3e50; border-bottom: 1px solid #eee; padding-bottom: 10px; margin-bottom: 15px; }}
            
            .grid-2 {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
            .grid-3 {{ display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 20px; }}
            
            /* Table Styling */
            table.dataTable thead th {{ background-color: #f1f3f5; color: #495057; }}
            .pass {{ color: #27ae60; font-weight: bold; background-color: #eafaf1; padding: 2px 8px; border-radius: 12px; font-size: 0.9em; }}
            .fail {{ color: #c0392b; font-weight: bold; background-color: #fdedec; padding: 2px 8px; border-radius: 12px; font-size: 0.9em; }}
            
            .metric-box {{ text-align: center; padding: 15px; background: #f8f9fa; border-radius: 8px; }}
            .metric-val {{ font-size: 24px; font-weight: bold; color: #2c3e50; }}
            .metric-label {{ font-size: 12px; color: #7f8c8d; text-transform: uppercase; letter-spacing: 1px; }}
        </style>
    </head>
    <body>
        <div class="header">
            <h1>🧪 AI-Driven Compound Analysis Report</h1>
            <p>Generated by Ahmed's Prediction Tool | {len(df)} Compounds Analyzed</p>
        </div>
        
        <div class="container">
            
            <!-- Summary Metrics -->
            <div class="grid-3">
                <div class="card metric-box">
                    <div class="metric-val">{len(df)}</div>
                    <div class="metric-label">Total Compounds</div>
                </div>
                <div class="card metric-box">
                    <div class="metric-val">{status_counts.get('PASS', 0)}</div>
                    <div class="metric-label">Passed bRo5</div>
                </div>
                <div class="card metric-box">
                    <div class="metric-val">{status_counts_lipinski.get('PASS', 0)}</div>
                    <div class="metric-label">Passed Lipinski</div>
                </div>
            </div>

            <!-- Top 5 Candidates -->
            <div class="card">
                <h2>🏆 Top 5 Candidates (Best Hybrid Score)</h2>
                <div style="overflow-x: auto;">
                    {df.head(5).to_html(classes='display compact', index=False, escape=False).replace('PASS', '<span class="pass">PASS</span>').replace('FAIL', '<span class="fail">FAIL</span>')}
                </div>
            </div>

            <!-- Charts Row 1 -->
            <div class="grid-2">
                <div class="card">
                    <h2>bRo5 Compliance</h2>
                    <div id="statusChart"></div>
                </div>
                <div class="card">
                    <h2>Lipinski Rule of 5 Compliance</h2>
                    <div id="lipinskiChart"></div>
                </div>
            </div>
            
            <!-- Charts Row 2 -->
            <div class="grid-2">
                <div class="card">
                    <h2>Chemical Space (MW vs LogP)</h2>
                    <div id="chemSpaceChart"></div>
                </div>
                <div class="card">
                    <h2>Similarity Distribution (Max Score per Compound)</h2>
                    <div id="simDistChart"></div>
                </div>
            </div>
            
            <!-- Top Docking Results (if applicable) -->
            {'<div class="card" style="margin-bottom: 20px;">' if has_docking else ''}
                {'<h2>🏆 Top 10 Docking Candidates</h2>' if has_docking else ''}
                {'<p style="color: #7f8c8d; margin-bottom: 15px;">Best binding affinity predictions (lower score = stronger binding)</p>' if has_docking else ''}
                {'<table style="width: 100%; border-collapse: collapse;">' if has_docking else ''}
                    {'<thead>' if has_docking else ''}
                        {'<tr style="background: #34495e; color: white;">' if has_docking else ''}
                            {'<th style="padding: 10px; text-align: left;">Rank</th>' if has_docking else ''}
                            {'<th style="padding: 10px; text-align: left;">Compound</th>' if has_docking else ''}
                            {'<th style="padding: 10px; text-align: left;">SMILES</th>' if has_docking else ''}
                            {'<th style="padding: 10px; text-align: right;">Hybrid Score</th>' if has_docking else ''}
                            {'<th style="padding: 10px; text-align: right;">Docking (kcal/mol)</th>' if has_docking else ''}
                            {f'<th style="padding: 10px; text-align: right;">Vina</th>' if has_docking and has_vina else ''}
                            {f'<th style="padding: 10px; text-align: right;">AutoDock-GPU</th>' if has_docking and has_autodock else ''}
                        {'</tr>' if has_docking else ''}
                    {'</thead>' if has_docking else ''}
                    {'<tbody>' if has_docking else ''}
                        {chr(10).join([f'''<tr style="border-bottom: 1px solid #ecf0f1;">
                            <td style="padding: 10px; font-weight: bold;">{i+1}</td>
                            <td style="padding: 10px;"><strong>{cand['name']}</strong></td>
                            <td style="padding: 10px; font-family: monospace; font-size: 11px;">{cand['smiles']}</td>
                            <td style="padding: 10px; text-align: right;">{cand['hybrid']}</td>
                            <td style="padding: 10px; text-align: right; color: #e74c3c; font-weight: bold;">{cand['docking']}</td>
                            {f'<td style="padding: 10px; text-align: right;">{cand.get("vina", "N/A")}</td>' if has_vina else ''}
                            {f'<td style="padding: 10px; text-align: right;">{cand.get("autodock", "N/A")}</td>' if has_autodock else ''}
                        </tr>''' for i, cand in enumerate(top_docking_candidates)]) if has_docking else ''}
                    {'</tbody>' if has_docking else ''}
                {'</table>' if has_docking else ''}
            {'</div>' if has_docking else ''}
            
            <!-- Charts Row 3 (Docking) -->
            <div class="grid-2" style="display: {'block' if has_docking else 'none'}; grid-template-columns: 1fr 1fr; gap: 20px; margin-bottom: 20px;">
                <div class="card" style="width: 100%;">
                     <h2>Docking Score Distribution (kcal/mol)</h2>
                     <div id="dockingHistChart"></div>
                </div>
                <div class="card" style="width: 100%;">
                     <h2>Consensus: Docking vs Hybrid Score</h2>
                     <div id="consensusChart"></div>
                </div>
            </div>

            <!-- Data Table -->
            <div class="card">
                <h2>Detailed Results Data</h2>
                <div style="overflow-x: auto;">
    """
    
    # Convert DataFrame to HTML (without index)
    # We remove the temporary Max columns before printing
    cols_to_drop = ['Max_Hybrid', 'Max_Tanimoto', 'Max_LLM']
    df_display = df.drop(columns=[c for c in cols_to_drop if c in df.columns])
    
    table_html = df_display.to_html(classes='display compact stripe hover', table_id='resultsTable', index=False, escape=False)
    
    # Highlight PASS/FAIL
    table_html = table_html.replace('PASS', '<span class="pass">PASS</span>')
    table_html = table_html.replace('FAIL', '<span class="fail">FAIL</span>')
    
    html_string += table_html
    
    html_string += f"""
                </div>
            </div>
        </div>
        
        <script>
            // Initialize DataTable
            $(document).ready( function () {{
                $('#resultsTable').DataTable({{
                    "pageLength": 10,
                    "scrollX": true,
                    "order": [] // Disable initial sort
                }});
            }});

            // 1. Pie Chart: bRo5 Status
            var statusData = {json.dumps(status_counts)};
            var pieData = [{{
                values: Object.values(statusData),
                labels: Object.keys(statusData),
                type: 'pie',
                marker: {{colors: ['#27ae60', '#c0392b', '#f39c12']}},
                hole: 0.4
            }}];
            Plotly.newPlot('statusChart', pieData, {{height: 350, margin: {{t:0, b:0, l:0, r:0}}}});

            // 1b. Pie Chart: Lipinski Status
            var statusDataLipinski = {json.dumps(status_counts_lipinski)};
            var pieDataLipinski = [{{
                values: Object.values(statusDataLipinski),
                labels: Object.keys(statusDataLipinski),
                type: 'pie',
                marker: {{colors: ['#27ae60', '#c0392b', '#f39c12']}},
                hole: 0.4
            }}];
            Plotly.newPlot('lipinskiChart', pieDataLipinski, {{height: 350, margin: {{t:0, b:0, l:0, r:0}}}});

            // 2. Scatter Plot: Chemical Space
            var rawData = {json.dumps(scatter_data)};
            var trace1 = {{
                x: rawData.map(d => d.x),
                y: rawData.map(d => d.y),
                mode: 'markers',
                type: 'scatter',
                text: rawData.map(d => d.name + '<br>' + d.status),
                marker: {{ 
                    size: 8,
                    color: rawData.map(d => d.status === 'PASS' ? '#27ae60' : '#c0392b'),
                    opacity: 0.6
                }}
            }};
            var layoutSpace = {{
                xaxis: {{title: 'Molecular Weight'}},
                yaxis: {{title: 'LogP'}},
                height: 350,
                margin: {{t:20, b:40, l:50, r:20}},
                hovermode: 'closest'
            }};
            Plotly.newPlot('chemSpaceChart', [trace1], layoutSpace);

            // 3. Histogram: Similarity Scores
            var simData = {json.dumps(sim_data)};
            if (simData.length > 0) {{
                var traceHybrid = {{
                    x: simData.map(d => d.hybrid),
                    type: 'histogram',
                    opacity: 0.6,
                    name: 'Hybrid Score',
                    marker: {{color: '#8e44ad'}}
                }};
                var traceTanimoto = {{
                    x: simData.map(d => d.tanimoto),
                    type: 'histogram',
                    opacity: 0.5,
                    name: 'Tanimoto',
                    marker: {{color: '#2980b9'}}
                }};
                var traceLLM = {{
                    x: simData.map(d => d.llm),
                    type: 'histogram',
                    opacity: 0.5,
                    name: 'LLM Cosine',
                    marker: {{color: '#e67e22'}}
                }};
                
                var layoutHist = {{
                    barmode: 'overlay',
                    xaxis: {{title: 'Similarity Score (0-1)'}},
                    yaxis: {{title: 'Count'}},
                    height: 350,
                    margin: {{t:20, b:40, l:50, r:20}}
                }};
                Plotly.newPlot('simDistChart', [traceHybrid, traceTanimoto, traceLLM], layoutHist);
            }}

            // 4. Docking Charts
            var hasDocking = { 'true' if has_docking else 'false' };
            if (hasDocking) {{
                var dockingRaw = {json.dumps(docking_data)};
                
                // Histogram
                var scores = dockingRaw.map(d => d.docking);
                var traceDock = {{
                    x: scores,
                    type: 'histogram',
                    marker: {{color: '#e74c3c'}},
                    name: 'Docking Score'
                }};
                var layoutDock = {{
                    xaxis: {{title: 'Binding Affinity (kcal/mol)'}},
                    yaxis: {{title: 'Frequency'}},
                    height: 350,
                    margin: {{t:20, b:40, l:50, r:20}}
                }};
                Plotly.newPlot('dockingHistChart', [traceDock], layoutDock);
                
                // Scatter (Consensus)
                var traceConsensus = {{
                    x: dockingRaw.map(d => d.hybrid),
                    y: dockingRaw.map(d => d.docking),
                    mode: 'markers',
                    type: 'scatter',
                    text: dockingRaw.map(d => d.name),
                    marker: {{
                        size: 9,
                        color: '#34495e',
                        opacity: 0.7
                    }}
                }};
                var layoutCon = {{
                    xaxis: {{title: 'Hybrid Similarity Score'}},
                    yaxis: {{title: 'Docking Score (kcal/mol)'}},
                    title: 'Lower Docking Score + Higher Hybrid Score = Better',
                    height: 350,
                    margin: {{t:40, b:40, l:50, r:20}},
                    hovermode: 'closest'
                }};
                Plotly.newPlot('consensusChart', [traceConsensus], layoutCon);
            }}
        </script>
    </body>
    </html>
    """
    
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_string)
    
    print(f"HTML Report generated at: {output_path}")
