import subprocess
import time
import json
import dash
import math
from dash import dcc, html
from dash.dependencies import Input, Output, State
import plotly.graph_objects as go
from network_layout import node_position

SOLVER_SCRIPT = "/home/nitishdumoliya/waterNetwork/wdnd/model/loop_wdn_heuristic.py"
JSON_DIR = "/home/nitishdumoliya/waterNetwork/wdnd/figure/json_file"
TOP_BAR_HEIGHT = "64px",
app = dash.Dash(__name__)
app.title = "Water Network Optimization Dashboard"
app.layout = html.Div(
    style={"width": "100%", "height": "100vh", "overflow": "hidden", "fontFamily": "Inter, Arial, sans-serif", "background": "#f3f4f6"},
    children=[
        html.H3("Water Distribution Network Optimization"),
        # html.Hr(style={"borderColor": "#e5e7eb"}),
        html.Div(
        style={"height": f"calc(100vh - {TOP_BAR_HEIGHT})",
            "display": "flex",
            "alignItems": "center",
            "justifyContent": "space-between",
            # "backgroundColor": "#e5e7eb",
            "padding": "8px",
            # "borderRadius": "10px",
            # "marginBottom": "10px"
        },
        children=[
            # ---------- LEFT SIDE (Controls) ----------
            html.Div(
                style={"display": "flex", "gap": "10px", "alignItems": "center"},
                children=[
                    html.Label("Select Network"),
                    dcc.Input(
                        id="data-number",
                        type="number",
                        min=1,
                        max=14,
                        step=1,
                        value=1,
                        style={"width": "70px"}
                    ),
                    html.Button(
                        "Run Optimization",
                        id="run-btn",
                        n_clicks=0
                    )
                ]
            ),
    
            # ---------- RIGHT SIDE (Output Info) ----------
            html.Div(
                id="top-info",
                style={
                    "fontWeight": "bold",
                    "whiteSpace": "nowrap",
                    "textAlign": "right"
                },
                children="Network: – | Objective: – | Time: –"
            )
        ]
    ),
        dcc.Loading(
            # type="circle",
            color="#2563eb",  # blue spinner
            children=dcc.Graph(
                id="wdn-graph",
                config={
                    "scrollZoom": True,
                    "doubleClick": "reset",
                    "displaylogo": False,
                    "modeBarButtonsToRemove": [
                        "lasso2d",
                        "select2d",
                        "autoScale2d",
                        "resetScale2d",
                        "toImage"
                    ]
                },
                style={
                    "height": "100vh",
                    "borderRadius": "12px"
                }
            )
        )
        # dcc.Loading(
        #     type="circle",
        #     children=dcc.Graph(id="wdn-graph",     
        #                        config={
        #                             "scrollZoom": True,        # ✅ mouse wheel zoom
        #                             # "displayModeBar": True,
        #                             "doubleClick": "reset", "displaylogo": False
        #                        }, 
        #                        style={"background": "white",
        #                                 "borderRadius": "10px",
        #                                 "padding": "10px",
        #                                 "boxShadow": "0 4px 10px rgba(0,0,0,0.15)"
        #                        }
        #     )
        # )
    ]
)

@app.callback(
    Output("wdn-graph", "figure"),
    Output("top-info", "children"),   # <- new output for top bar
    Input("run-btn", "n_clicks"),
    State("data-number", "value"),
    prevent_initial_call=True
)

def run_solver(n_clicks, data_number):

    # Run solver as separate process
    subprocess.run(
        ["python3", SOLVER_SCRIPT, str(data_number)],
        check=True
    )

    # Load JSON after solver finishes
    json_file = f"{JSON_DIR}/solution_{data_number-1}.json"
    node_pos = node_position[data_number-1]

    return build_plot(json_file, node_pos, data_number-1)

def build_plot(solution_file, node_pos, data_number):
    # --------------------------------------------------
    # Load solution
    # --------------------------------------------------
    with open(solution_file) as f:
        sol = json.load(f)
    def parse_key(k):
        return tuple(map(int, k.replace("(", "").replace(")", "").split(",")))
    q = {parse_key(k): v for k, v in sol["q"].items()}
    L = {parse_key(k): v for k, v in sol.get("L", {}).items()}
    # Pipe info l[i,j,d]
    pipe_info = {}
    for k_str, val in sol.get("l", {}).items():
        if val < 1e-6:
            continue
        i, j, d = map(int, k_str.replace("(", "").replace(")", "").split(","))
        pipe_info.setdefault((i, j), []).append((d, val))
    h    = {int(k): v for k, v in sol["h"].items()}
    hmin = {int(k): v for k, v in sol["hmin"].items()}
    D    = {int(k): v for k, v in sol["D"].items()}
    # source = sol["source"]
    source = sol.get("source", [])
    # --------------------------------------------------
    # Figure scaling + node radius
    # --------------------------------------------------
    node_marker_size = 20
    BASE = 1800
    xs = [p[0] for p in node_pos.values()]
    ys = [p[1] for p in node_pos.values()]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    data_w = xmax - xmin
    data_h = ymax - ymin
    if data_w >= data_h:
        fig_w = BASE
        fig_h = int(BASE * data_h / data_w)
    else:
        fig_h = BASE
        fig_w = int(BASE * data_w / data_h)
    node_radius = (node_marker_size / 2) * (data_w / fig_w)
    # --------------------------------------------------
    # Diameter → thickness + color (DISCRETE)
    # --------------------------------------------------
    arc_diameter = {}
    # all_diams = sorted({
    #     max(d for d, _ in pipes)
    #     for pipes in pipe_info.values()
    # })
    all_diams = sorted({ d for pipes in pipe_info.values() for d, _ in pipes })
    # Visually separated thickness levels
    thickness_levels = [3, 6, 9, 12, 15]
    # Distinct categorical colors
    color_palette = [
        "#1f77b4", "#ff7f0e", "#2ca02c",
        "#d62728", "#9467bd", "#8c564b",
        "#e377c2", "#7f7f7f", "#bcbd22",
        "#17becf"
    ]
    diam_to_style = {}
    for i, d in enumerate(all_diams):
        diam_to_style[d] = dict(
            width=thickness_levels[min(i, len(thickness_levels)-1)],
            color=color_palette[i % len(color_palette)]
        )
    for (i, j), pipes in pipe_info.items():
        arc_diameter[(i, j)] = max(d for d, _ in pipes)
    # --------------------------------------------------
    # Containers
    # --------------------------------------------------
    edge_groups = {}  # d → lines
    arrows = []
    click_x, click_y, click_text = [], [], []
    flow_tx, flow_ty, flow_txt = [], [], []
    pipe_tx, pipe_ty, pipe_txt = [], [], []
    # --------------------------------------------------
    # EDGES (arc hover + multi-diameter coloring)
    # --------------------------------------------------
    for (i, j), flow in q.items():
    
        if i not in node_pos or j not in node_pos:
            continue
        x0, y0 = node_pos[i]
        x1, y1 = node_pos[j]
        xs, ys, xe, ye = (x0, y0, x1, y1) if flow >= 0 else (x1, y1, x0, y0)
        dx, dy = xe - xs, ye - ys
        dist = math.hypot(dx, dy)
        if dist < 1e-9:
            continue
        sx = xs + node_radius * dx / dist
        sy = ys + node_radius * dy / dist
        ex = xe - node_radius * dx / dist
        ey = ye - node_radius * dy / dist
        mx, my = (sx + ex) / 2, (sy + ey) / 2
        # --------------------------------------------------
        # Hover info (UNCHANGED – arc level)
        # --------------------------------------------------
        pipes = pipe_info.get((i, j), [])
        pipe_label = ", ".join([f"D{d}:{l:.1f}" for d, l in pipes])
        hover = (
            f"<b>Arc {i} → {j}</b><br>"
            f"Flow: {flow:.5f}<br>"
            f"Length: {L.get((i,j),0):.2f}<br>"
            f"<b>Pipes:</b><br>{pipe_label}"
        )
        click_x.append(mx)
        click_y.append(my)
        click_text.append(hover)
        flow_tx.append(mx)
        flow_ty.append(my)
        flow_txt.append(f"{flow:.3f}")
        pipe_tx.append(mx)
        pipe_ty.append(my - 0.02 * data_h)
        pipe_txt.append(pipe_label)
        # --------------------------------------------------
        # Multi-diameter segmentation
        # --------------------------------------------------
        total_len = sum(l for _, l in pipes)
        if total_len <= 0:
            continue
        cur_len = 0.0
        for d, l in pipes:
            frac_s = cur_len / total_len
            frac_e = (cur_len + l) / total_len
            seg_sx = sx + frac_s * dx
            seg_sy = sy + frac_s * dy
            seg_ex = sx + frac_e * dx
            seg_ey = sy + frac_e * dy
            style = diam_to_style[d]
            edge_groups.setdefault(d, {"x": [], "y": []})
            edge_groups[d]["x"] += [seg_sx, seg_ex, None]
            edge_groups[d]["y"] += [seg_sy, seg_ey, None]
            cur_len += l
        # --------------------------------------------------
        # Arrow (single arrow for whole arc)
        # --------------------------------------------------
        d_arrow = max(d for d, _ in pipes)
        style = diam_to_style[d_arrow]
        arrows.append(dict(
            ax=sx, ay=sy, x=ex, y=ey,
            xref="x", yref="y",
            axref="x", ayref="y",
            showarrow=True,
            arrowhead=3,
            arrowsize=2,
            arrowwidth=1,
            arrowcolor= "black"
            # arrowcolor= style["color"]
        ))
    # --------------------------------------------------
    # Traces
    # --------------------------------------------------
    traces = []
    for d in sorted(edge_groups.keys()):   # ← increasing order
        data = edge_groups[d]
        style = diam_to_style[d]
        traces.append(go.Scatter(
            x=data["x"], y=data["y"],
            mode="lines",
            line=dict(
                width=style["width"],
                color=style["color"]
            ),
            name=f"Diameter D{d}",
            hoverinfo="skip"
        )) 
    traces.append(go.Scatter(
        x=click_x, y=click_y,
        mode="markers",
        marker=dict(size=0, opacity=0),
        hoverinfo="text",
        hovertext=click_text,
        showlegend=False
    ))
    # --------------------------------------------------
    # Nodes
    # --------------------------------------------------
    nx, ny, ntext, nlabel, ncolor = [], [], [], [], []
    for n, (x, y) in node_pos.items():
        nx.append(x); ny.append(y)
        nlabel.append(str(n))
        ntext.append(
            f"<b>Node {n}</b><br>"
            f"Head: {h.get(n,0):.2f}<br>"
            f"Min Head: {hmin.get(n,0):.2f}<br>"
            f"Demand: {D.get(n,0):.5f}"
        )
        ncolor.append("royalblue" if n in source else "skyblue")
    traces.append(go.Scatter(
        x=nx, y=ny,
        mode="markers+text",
        text=nlabel,
        textposition="middle center",
        hoverinfo="text",
        hovertext=ntext,
        marker=dict(size=node_marker_size, color=ncolor, line=dict(width=1.5, color="black")),
        name="Nodes"
    ))
    # --------------------------------------------------
    # Figure
    # --------------------------------------------------
    data_list = [
        "d1_bessa",
        "d2_shamir",
        "d3_hanoi",
        "d4_double_hanoi",
        "d5_triple_hanoi",
        "d6_newyork",
        "d7_blacksburg",
        "d8_fossolo_iron",
        "d9_fossolo_poly_0",
        "d10_fossolo_poly_1",
        "d11_kadu",
        "d12_pescara",
        "d13_modena",
        "d14_balerma"
    ]
    fig = go.Figure(traces)
    # fig.update_layout(
    #     title=f"WDN = {data_list[data_number]}, Total Cost = {sol['objective']:.2f}",
    #     annotations=arrows,
    #     dragmode="pan",
    #     hovermode="closest",
    #     autosize=False,
    #     xaxis=dict(fixedrange=False),
    #     yaxis=dict(scaleanchor="x", fixedrange=False),
    #     plot_bgcolor="white",
    #     paper_bgcolor="white",
    #     width=fig_w,
    #     height=fig_h,
    #     # margin=dict(l=20, r=20, t=50, b=20),
    # )
    fig.update_layout(
        # title=f"WDN = {data_list[data_number]}, Total Cost = {sol['objective']:.2f}",
        annotations=arrows,
        dragmode="pan",
        hovermode="closest",
        # hovermode="x unified",
        autosize=False,
        width=1900,
        height=950,
        xaxis=dict(fixedrange=False),
        yaxis=dict(scaleanchor="x", fixedrange=False),
        plot_bgcolor="white",
        paper_bgcolor="white",
    )
    fig.write_html(
        f"../figure/json_file/wdn_interactive_solution{data_number+1}.html",
        auto_open=False,
        config={"scrollZoom": True, 
                }
    )
    network_info = f"Network: {data_list[data_number]} | Objective: {sol['objective']:.2f} | Time: {sol['solve_time']:.2f} s"
    print("✔ Diameter-colored & thickness-separated WDN plot saved")
    return fig, network_info

if __name__ == "__main__":
    app.run(debug=True)
