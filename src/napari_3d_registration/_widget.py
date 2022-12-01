"""
This module is an example of a barebones QWidget plugin for napari

It implements the Widget specification.
see: https://napari.org/stable/plugins/guides.html?#widgets

Replace code below according to your needs.
"""
from typing import TYPE_CHECKING

from magicgui import widgets
from qtpy.QtWidgets import QVBoxLayout, QPushButton, QWidget, QTabWidget
from registrationtools import SpatialRegistration, TimeRegistration
from pathlib import Path
from IO import imread

class TimeRegistrationWidget(QWidget):
    @property
    def parameters(self):
        self._parameters = {
            "path_to_data": self.path_data_file.value,
            "file_name": self.file_format.value,
            "trsf_folder": self.trsf_folder.value,
            "output_format": self.output_path.value,
            "projection_path": self.projection_path.value,
            "check_TP": 0,
            "path_to_bin": "",
            "voxel_size": [
                self.x_vox.value,
                self.y_vox.value,
                self.z_vox.value,
            ],
            "first": self.first_tp.value,
            "last": self.last_tp.value,
            "not_to_do": self.no_tp.value,
            "compute_trsf": self.compute_trsf.value,
            "ref_TP": self.ref_tp.value,
            "trsf_type": self.trsf_type.value,
            "padding": self.padding.value,
            "recompute": self.recompute.value,
            "pre_2D": self.pre2D.value,
            "low_th": self.low_th.value,
            "registration_depth_end": self.reg_depth_end,
            "registration_depth_start": self.reg_depth_start,
            "apply_trsf": self.apply_trsfs.value,
            "time_tag": "",
            "out_bdv": "",
        }
        return self._parameters

    def make_file_search(self, label, filter, where=None, default=None):
        label = widgets.Label(value=label)
        if where is None:
            where = label
        if default is None:
            default = Path(".").absolute()
        box = widgets.FileEdit(value=default, filter=filter)
        if isinstance(where, str):
            self.__dict__[where] = box
        else:
            where = box
        container = widgets.Container(widgets=[label, box], labels=False)
        return container

    def make_text_edit(self, label, where=None, default=None):
        label = widgets.Label(value=label)
        if where is None:
            where = label
        if default is None:
            default = Path(".").absolute()
        box = widgets.LineEdit(value=default)
        if isinstance(where, str):
            self.__dict__[where] = box
        else:
            where = box
        container = widgets.Container(widgets=[label, box], labels=False)
        return container

    def make_tick_box(self, label, default=True, where=None):
        label = widgets.Label(value=label)
        if where is None:
            where = label
        tick = widgets.CheckBox(value=default)
        if isinstance(where, str):
            self.__dict__[where] = tick
        else:
            where = tick
        container = widgets.Container(
            widgets=[label, tick], labels=False, layout="horizontal"
        )
        return container

    def make_button(self, name, target):
        button = QPushButton(name)
        button.clicked.connect(target)
        button.native = button
        return button

    def __init__(self, napari_viewer):
        super().__init__()
        self.viewer = napari_viewer

        ## json file parameterization:
        json = self.make_file_search("Parameter file", "*.json", "json_file")
        start_reg = self.make_button("Run!", self._on_click_file)
        file_params = widgets.Container(
            widgets=[json, start_reg], labels=False
        )
        file_params.native.layout().addStretch(1)


        ## Manual parameterization
        ### Paths definition
        path_data = self.make_file_search(
            label="Path to images", filter=None, where="path_data_file"
        )
        file_name = self.make_text_edit(
            label="File format",
            where="file_format",
            default="image{t:03d}.tiff",
        )
        trsf_folder = self.make_file_search(
            label="Trsf folder", filter="", where="trsf_folder"
        )
        output_path = self.make_file_search(
            label="Output Format",
            filter="",
            where="output_path",
            default="_registered",
        )
        projection_path = self.make_file_search(
            label="Path to projection", filter="", where="projection_path"
        )
        path_tab = widgets.Container(
            widgets=[
                path_data,
                file_name,
                trsf_folder,
                output_path,
                projection_path
            ],
            labels=False
        )
        ### Geometry definition
        vox_label = widgets.Label(value="Voxel size")
        x_label = widgets.Label(value='x')
        self.x_vox = widgets.FloatText(value=1)
        y_label = widgets.Label(value='y')
        self.y_vox = widgets.FloatText(value=1)
        z_label = widgets.Label(value='z')
        self.z_vox = widgets.FloatText(value=6)
        resolution = widgets.Container(
            widgets=[
                x_label,
                self.x_vox,
                y_label,
                self.y_vox,
                z_label,
                self.z_vox,
            ],
            layout="vertical",
            labels=False,
        )
        vox_res = widgets.Container(
            widgets=[vox_label, resolution], labels=False
        )
        geometry_tab = widgets.Container(
            widgets=[
                vox_res
            ],
            labels=False
        )

        ### Time definition
        first_tp_label = widgets.Label(value="First time point:")
        self.first_tp = widgets.IntText(value=0)
        f_tp = widgets.Container(
            widgets=[first_tp_label, self.first_tp],
            layout="horizontal",
            labels=False,
        )

        last_tp_label = widgets.Label(value="Last time point:")
        self.last_tp = widgets.IntText(value=100)
        l_tp = widgets.Container(
            widgets=[last_tp_label, self.last_tp],
            layout="horizontal",
            labels=False,
        )

        ref_tp_label = widgets.Label(value="Reference time point:")
        self.ref_tp = widgets.IntText(value=50)
        r_tp = widgets.Container(
            widgets=[ref_tp_label, self.ref_tp],
            layout="horizontal",
            labels=False,
        )

        time_tab = widgets.Container(
            widgets=[
                f_tp,
                l_tp,
                r_tp
            ],
            labels=False
        )
        
        ### Trsf parameterization
        trsf_type_label = widgets.Label(value="Transformation type")
        self.trsf_type = widgets.ComboBox(
            value="rigid", choices=["rigid", "affine"]
        )
        trsf_type = widgets.Container(
            widgets=[trsf_type_label, self.trsf_type], labels=False
        )
        
        compute_trsf = self.make_tick_box(
            "Compute trsf", True, where="compute_trsf"
        )
        padding = self.make_tick_box("Padding", True, where="padding")
        recompute = self.make_tick_box("Force recompute", True, where="recompute")
        apply_trsf = self.make_tick_box(
            "Apply Trsf", True, where="apply_trsfs"
        )
        
        trsf_tab = widgets.Container(
            widgets=[
                trsf_type,
                compute_trsf,
                padding,
                apply_trsf,
                recompute
            ],
            labels=False
        )

        ### Advanced parameterization
        no_tp = self.make_text_edit("Time points to skip:", "no_tp", "[]")
        pre2D = self.make_tick_box("Pre 2D registration", True, where="pre2D")
        low_th = self.make_tick_box("Low threshold", True, where="low_th")
        low_th_value_label = widgets.Label(value="Low threshold value:")
        self.low_th_val = widgets.IntText(value=250)
        low_th_value = widgets.Container(
            widgets=[low_th_value_label, self.low_th_val],
            layout="horizontal",
            labels=False,
        )
        reg_depth_label_start = widgets.Label(value="Registration depth start:")
        self.reg_depth_start = widgets.IntText(value=3)
        r_depth_start = widgets.Container(
            widgets=[reg_depth_label_start, self.reg_depth_start],
            layout="horizontal",
            labels=False,
        )
        reg_depth_label_end = widgets.Label(value="Registration depth end:")
        self.reg_depth_end = widgets.IntText(value=3)
        r_depth_end = widgets.Container(
            widgets=[reg_depth_label_end, self.reg_depth_end],
            layout="horizontal",
            labels=False,
        )

        advanced_tab = widgets.Container(
            widgets=[
                no_tp,
                pre2D,
                low_th,
                low_th_value,
                r_depth_start,
                r_depth_end,
            ],
            labels=False
        )
        
        plot_trsf = self.make_tick_box("Plot trsf before apply", True, where="plot_trsf")
        start_reg_manual = self.make_button("Run!", self._on_click_manual)

        tab_controls = QTabWidget()
        tab_controls.addTab(file_params.native, "From File")

        sub_tab = QTabWidget()
        sub_tab.addTab(path_tab.native, "Paths")
        sub_tab.addTab(trsf_tab.native, "Trsf")
        sub_tab.addTab(geometry_tab.native, "Geometry")
        sub_tab.addTab(time_tab.native, "Time")
        sub_tab.addTab(advanced_tab.native, "Advanced")
        sub_tab.native = sub_tab
        manual_controler = widgets.Container(
            widgets=[sub_tab, plot_trsf, start_reg_manual], labels=False
        )
        # sub_tab.addStretch(1)
        tab_controls.addTab(manual_controler.native, "Manual")
        # tab_controls.addStretch(1)

        layout = QVBoxLayout()
        layout.addStretch(1)
        self.setLayout(layout)
        self.layout().addWidget(tab_controls)

    def _on_click_manual(self):
        tr = TimeRegistration(self.parameters)
        tr.run_trsf()
        for p in tr.params:
            if p.projection_path is not None:
                p_to_data = p.projection_path
                num_s = p.file_name.find("{")
                num_e = p.file_name.find("}") + 1
                f_name = p.file_name.replace(p.file_name[num_s:num_e], "")
                im = imread(
                    p_to_data + f_name.replace(p.im_ext, "xyProjection.tif")
                ).transpose(2, 1, 0)
                self.viewer.add_image(im)

    def _on_click_file(self):
        tr = TimeRegistration(self.json_file.value)
        tr.run_trsf()
        p = tr.params[0]
        if p.projection_path is not None:
            p_to_data = p.projection_path
            num_s = p.file_name.find("{")
            num_e = p.file_name.find("}") + 1
            f_name = p.file_name.replace(p.file_name[num_s:num_e], "")
            im = imread(
                p_to_data + f_name.replace(p.im_ext, "xyProjection.tif")
            ).transpose(2, 1, 0)
            self.viewer.add_image(im)
