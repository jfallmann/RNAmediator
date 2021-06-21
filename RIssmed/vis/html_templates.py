import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from typing import List, Union
import os
import base64
import zipfile







def search_inputs():
    columns = "p-1 col-md-3 col-10 d-flex justify-content-center align-items-center text-center"
    menu = html.Div(
        [
            html.Div(
                html.Button("", className="prev-button", id="prev-button", n_clicks=0),
                className="col-2 col-md-1 text-center order-2 order-md-1 m-1 m-md-0"),
            html.Div(html.Div([
                html.Nobr([
                    html.Label(
                        "Chr",
                        htmlFor="search-chr-input",
                        className="search-label"

                    ),
                    dcc.Input(id="search-chr-input",
                              className="search-input",
                              placeholder="Search Chromosome",
                              )
                ], className=columns),
                html.Nobr([
                    html.Label(
                        "GOI",
                        htmlFor="search-goi-input",
                        className="search-label"
                    ),
                    dcc.Input(id="search-goi-input",
                              placeholder="Search Gene",
                              className="search-input",
                              )
                ], className=columns),
                html.Nobr([
                    html.Label(
                        "Start",
                        htmlFor="search-span-start-input",
                        className="search-label"
                    ),
                    dcc.Input(id="search-span-start-input",
                              placeholder="Start",
                              className="search-input",
                              )
                ], className=columns),
                html.Nobr([
                    html.Label(
                        "End",
                        htmlFor="search-span-end-input",
                        className="search-label"
                    ),
                    dcc.Input(id="search-span-end-input",
                              placeholder="End",
                              className="search-input",
                              )
                ], className=columns),
            ], className="row justify-content-center justify-content-md-between"),
                className="col-12 col-md-10 order-md-2 order-1"),
            html.Div(html.Button("", className="next-button", id="next-button", n_clicks=0),
                     className="col-2 col-md-1 text-center order-3 order-md-3 m-1 m-md-0")
        ],
        className="row justify-content-around align-items-center p-2",

    )
    return menu


def interesting_table(interesting, prev_clicks: int = 0, next_clicks: int = 0,
                      sorting: str = "Max_Value", sorting_clicks: int = 0, number_of_interesting: int = 10):
    page = next_clicks - prev_clicks
    if len(interesting) == 0:
        header = ["WARNING", "No", "matching", "entries"]
        interesting = [(0, "-", "-", "-", "-")]
        clickable = False
    else:
        clickable = True
        header = interesting[0].keys()
        interesting = [(x + 1 + page * number_of_interesting, *element) for x, element in enumerate(interesting)]
    menu = [html.Div([
                      html.Div(html.Table(
                          [html.Tr([html.Td(html.Button("#", style={"width": "100%", "pointer-events": "None"}),
                                            className="interesting-table-header")] +
                                   list(table_header_generator(header, sort=sorting, sorting_clicks=sorting_clicks,
                                                               clickable=clickable))
                                   )] +
                          [html.Tr(
                              list(tablerow_generator(entry))
                          ) for y, entry in enumerate(interesting)],
                          className="interesting-table-tablerows",
                          style={"width": "100%", "margin": "auto", 'text-align': 'center'}

                      ), style={"width": "100%", "overflow": "auto"}, className="interesting-table col-12 m-0 p-0"),
                      ],
                     className="col-11 justify-content-center")]
    return menu



def table_header_generator(row, sort, sorting_clicks, clickable: bool = True):
    for element in row:
        if clickable:
            style = {"white-space": "nowrap", "margin": "auto", "width": "100%",
                     "font-weight": "bold"}
        else:
            style = {"white-space": "nowrap", "margin": "auto", "width": "100%",
                     "font-weight": "bold", "pointer-events": "None", "color": "#FF7676"}
        if element == sort:
            if sorting_clicks % 2:
                order = "\U000025B2"
            else:
                order = "\U000025BC"

        else:
            order = ""
        column = html.Td(html.Button(f"{element.replace('_', ' ')}{order}", n_clicks=sorting_clicks,
                                     id={"index": "sorting",
                                         "type": "interesting-table-header-button",
                                         "name": f"{element}"},
                                     style=style),
                         className="interesting-table-header",
                         style={"font-weight": "bold"})
        yield column


def tablerow_generator(row):
    for x, col in enumerate(row):
        if x == 0:
            style = {"width": "7%"}
            classname = "interesting-table-index-column"
        elif x == 1:
            style = {"width": "10%"}
            classname = "interesting-table-column"
        else:
            style = {"width": f"{(100 - 17) / (len(row) - 2)}%"}
            classname = "interesting-table-column"
        if type(col) == float:
            col_to_show = round(col, 3)
        else:
            col_to_show = col
            col = ""
        column = html.Td(html.Button(col_to_show, id={"index": f"{x}-{row[1]}-{row[2]}",
                                                      "type": "interesting-table-button",
                                                      "chrom": f"{row[1]}", "goi": row[2],
                                                      "start": row[3],
                                                      "end": row[4]},
                                     n_clicks=0,
                                     className="interesting-table-button",
                                     title=f"{col}",
                                     style={"white-space": "nowrap", "margin": "auto", "width": "100%"}),

                         className=classname, style=style)
        yield column


def modal_image_download():
    modal = dbc.Modal([
                              dbc.ModalHeader("Select Format"),
                              dbc.ModalBody([
                                  html.Div([
                                      html.Button("SVG", id="svg-download", className="btn btn-primary col-3"),
                                      html.Button("PNG", id="png-download", className="btn btn-primary col-3")
                                  ], className="row justify-content-around")
                              ]),
                              dbc.ModalFooter(
                                  dbc.Button(
                                      "Close", id="close", className="ml-auto", n_clicks=0
                                  )
                              ),
                          ], id="modal")
    return modal


def get_ingo(first_name, name, assets_dir):
    if first_name in ["ingo", "Ingo"] and name in ["Flamingo", "flamingo"]:
        img_path = os.path.join(assets_dir, "animation.zip")
        with zipfile.ZipFile(img_path) as file:
            img = file.read("animation_2.svg", pwd=b"ingo")
        encoded_img = base64.b64encode(img)
        svg = 'data:image/svg+xml;base64,{}'.format(encoded_img.decode())
        ingo = html.Div(html.Div(
            [
                html.H3("Say Hello to Ingo", style={"text-align": "center"}),
                html.Div(html.Img(src=svg, style={"margin": "auto"}), style={"margin": "auto", "margin-top": "20px",
                                                                             "text-align": "center"})
            ], className="databox"), className="col-12 p-1 justify-content-center", id="ingo-box")
    else:
        ingo = []
    return ingo