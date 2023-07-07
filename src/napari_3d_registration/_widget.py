"""
This module is an example of a barebones QWidget plugin for napari

It implements the Widget specification.
see: https://napari.org/stable/plugins/guides.html?#widgets

Replace code below according to your needs.
"""
from abc import abstractproperty, abstractmethod
from typing import TYPE_CHECKING

from magicgui import widgets
from qtpy.QtWidgets import (
    QVBoxLayout,
    QPushButton,
    QWidget,
    QTabWidget,
    QFileDialog,
    QComboBox,
    QStackedWidget,
    QMessageBox,
)
from registrationtools import SpatialRegistration, TimeRegistration
from pathlib import Path
from IO import imread
from shutil import which
import json


class RegistrationWidget(QWidget):
    @abstractproperty
    def link_to_parameters(self):
        self._link_to_parameters = {}
        return self._link_to_parameters

    @property
    def parameters(self):
        self._parameters = {}
        for k, v in self.link_to_parameters.items():
            ## Case of list of parameters that have a value attribute
            if isinstance(v, list) and all([hasattr(vs, "value") for vs in v]):
                val = []
                for vs in v:
                    val.append(vs.value)
            else:
                if hasattr(v, "value"):
                    val = v.value
                else:
                    val = v

                if isinstance(val, Path):
                    val = str(val)
            self._parameters[k] = val
        return self._parameters

    def make_file_search(
        self, label, filter, where=None, default=None, mode="r"
    ):
        label = widgets.Label(value=label)
        if where is None:
            where = label
        if default is None:
            default = Path(".").absolute()
        box = widgets.FileEdit(value=default, filter=filter, mode=mode)
        if isinstance(where, str):
            self.__dict__[where] = box
        else:
            where = box
        container = widgets.Container(widgets=[label, box], labels=False)
        return container

    def make_text_edit(self, label, where=None, default=None, eval_val=False):
        label = widgets.Label(value=label)
        if where is None:
            where = label
        if default is None:
            default = Path(".").absolute()
        if eval_val:
            box = widgets.LiteralEvalLineEdit(value=default)
        else:
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
        button.name = name
        return button

    def make_voxel_label(self, label="Voxel size", where="{axis}_vox"):
        vox_label = widgets.Label(value=label)
        x_label = widgets.Label(value="x")
        self.__dict__[where.format(axis="x")] = widgets.FloatText(
            value=1, min=0, max=100, step=1e-5
        )
        cont_x = widgets.Container(
            widgets=[x_label, self.__dict__[where.format(axis="x")]],
            layout="horizontal",
            labels=False,
        )
        y_label = widgets.Label(value="y")
        self.__dict__[where.format(axis="y")] = widgets.FloatText(
            value=1, min=0, max=100, step=1e-5
        )
        cont_y = widgets.Container(
            widgets=[y_label, self.__dict__[where.format(axis="y")]],
            layout="horizontal",
            labels=False,
        )
        z_label = widgets.Label(value="z")
        self.__dict__[where.format(axis="z")] = widgets.FloatText(
            value=6, min=0, max=100, step=1e-5
        )
        cont_z = widgets.Container(
            widgets=[z_label, self.__dict__[where.format(axis="z")]],
            layout="horizontal",
            labels=False,
        )
        resolution = widgets.Container(
            widgets=[
                cont_x,
                cont_y,
                cont_z,
            ],
            layout="vertical",
            labels=False,
        )
        vox_res = widgets.Container(
            widgets=[vox_label, resolution], labels=False
        )
        return vox_res

    def _propagate_json(self):
        with open(self.json_file.value) as f:
            json_params = json.load(f)
        if "voxel_size" in json_params and not "voxel_size_out" in json_params:
            json_params["voxel_size_out"] = json_params["voxel_size"]
        if "ref_voxel" in json_params and not "out_voxel" in json_params:
            json_params["out_voxel"] = json_params["ref_voxel"]
        for k, v in json_params.items():
            ## Classic parameters
            if k in self.link_to_parameters and hasattr(
                self.link_to_parameters[k], "value"
            ):
                if isinstance(v, list):
                    v = v[0]
                self.link_to_parameters[k].value = v
            ## Parameters that are lists of parameters
            elif isinstance(self.link_to_parameters[k], list) and all(
                [hasattr(vi, "value") for vi in self.link_to_parameters[k]]
            ):
                for i, vs in enumerate(v):
                    self.link_to_parameters[k][i].value = vs
            ## More complex parameters with dedicated getters and setters
            else:
                self.__setattr__(k, v)
                self.link_to_parameters[k] = self.__getattribute__(k)

    def _save_json(self):
        options = QFileDialog.Options()
        fileName, _ = QFileDialog.getSaveFileName(
            self,
            "Save parameter file",
            "parameters",
            "json Files (*.json)",
            options=options,
        )
        if fileName:
            with open(fileName, "w") as f:
                tmp = self.parameters.copy()
                if "trsf_paths" in tmp:
                    tmp["trsf_paths"] = [tmp["trsf_paths"]] * len(
                        tmp["flo_ims"]
                    )
                json.dump(tmp, f, indent=2)

    @abstractmethod
    def make_manual_parameterization(self):
        manual_controler = widgets.Label(value="Hello World")
        return manual_controler

    @abstractmethod
    def _on_click_manual(self):
        print("Manual Click")

    @abstractmethod
    def _on_click_file(self):
        print("File Click")

    @abstractmethod
    def _check_binaries(self):
        return []

    def __init__(self, napari_viewer):
        super().__init__()
        self.viewer = napari_viewer

        binaries_not_found = self._check_binaries()
        if 0 < len(binaries_not_found):
            not_found = "\n\t".join(binaries_not_found)
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Missing binaries")
            msg.setInformativeText(
                (
                    "The following binaries necessary to run this plugin were not found:\n"
                    "\t" + not_found
                )
            )
            msg.setWindowTitle("Missing binaries")
            msg.exec_()
            return

        ## json file parameterization:
        json = self.make_file_search("Parameter file", "*.json", "json_file")
        json.changed.connect(self._propagate_json)
        clear_tick = self.make_tick_box("Clear\nImages", default=False, where="clear_file")
        start_reg = self.make_button("Run!", self._on_click_file)
        run_box = widgets.Container(
            widgets=[clear_tick, start_reg],
            labels=False,
            layout="horizontal"
        )
        self.file_params = widgets.Container(
            widgets=[json, run_box], labels=False
        )
        self.file_params.native.layout().addStretch(1)

        ## Manual parameterization
        manual_controler = self.make_manual_parameterization()

        self.tab_controls = QTabWidget()
        self.tab_controls.addTab(self.file_params.native, "From File")
        self.tab_controls.addTab(manual_controler.native, "Manual")

        layout = QVBoxLayout()
        layout.addStretch(1)
        layout.setSpacing(0)
        self.setLayout(layout)
        self.layout().addWidget(self.tab_controls)


class TimeRegistrationWidget(RegistrationWidget):
    @property
    def link_to_parameters(self):
        self._link_to_parameters = {
            "path_to_data": self.path_data_file,
            "file_name": self.file_format,
            "trsf_folder": self.trsf_folder,
            "output_format": self.output_path,
            "projection_path": self.projection_path,
            "check_TP": False,
            "path_to_bin": "",
            "voxel_size": [
                self.x_vox,
                self.y_vox,
                self.z_vox,
            ],
            "voxel_size_out": [
                self.x_vox_out,
                self.y_vox_out,
                self.z_vox_out,
            ],
            "first": self.first_tp,
            "last": self.last_tp,
            "not_to_do": self.no_tp,
            "compute_trsf": self.compute_trsf,
            "ref_TP": self.ref_tp,
            "trsf_type": self.trsf_type,
            "padding": self.padding,
            "recompute": self.recompute,
            "pre_2D": self.pre2D,
            "low_th": self.low_th,
            "registration_depth_end": self.reg_depth_end,
            "registration_depth_start": self.reg_depth_start,
            "apply_trsf": self.apply_trsfs,
            "time_tag": "",
            "out_bdv": "",
            "plot_trsf": self.plot_trsf,
        }
        return self._link_to_parameters

    @property
    def low_th(self):
        if not self.low_th_bool.value:
            self._low_th = 0
        else:
            self._low_th = self.low_th_val.value
        return self._low_th

    @low_th.setter
    def low_th(self, value):
        self._low_th = value
        self.low_th_val.value = int(value if value else 0)
        self.low_th_bool.value = bool(value != 0)

    def _on_click_manual(self):
        if self.clear_manual:
            for l in self.viewer.layers[:]:
                self.viewer.layers.remove(l)
        tr = TimeRegistration(self.parameters)
        tr.run_trsf()
        for p in tr.params:
            if p.projection_path is not None:
                p_to_data = p.projection_path
                num_s = p.file_name.find("{")
                num_e = p.file_name.find("}") + 1
                f_name = p.file_name.replace(p.file_name[num_s:num_e], "")
                proj_name = f_name.replace(p.im_ext, "{:s}Projection.tif")
                im_xy = imread(p_to_data + proj_name.format("xy")).transpose(
                    2, 1, 0
                )
                im_xz = imread(p_to_data + proj_name.format("xz")).transpose(
                    2, 1, 0
                )
                im_yz = imread(p_to_data + proj_name.format("yz")).transpose(
                    2, 1, 0
                )
                self.viewer.add_image(
                    im_xy,
                    name=proj_name.format("xy"),
                    scale=(1, p.voxel_size_out[1], p.voxel_size_out[0]),
                )
                self.viewer.add_image(
                    im_yz,
                    name=proj_name.format("xz"),
                    scale=(1, p.voxel_size_out[2], p.voxel_size_out[1]),
                    visible=False,
                )
                self.viewer.add_image(
                    im_xz,
                    name=proj_name.format("yz"),
                    scale=(1, p.voxel_size_out[2], p.voxel_size_out[0]),
                    visible=False,
                )

    def _on_click_file(self):
        if self.clear_file:
            for l in self.viewer.layers[:]:
                self.viewer.layers.remove(l)
        tr = TimeRegistration(self.json_file.value)
        tr.run_trsf()
        p = tr.params[0]
        if p.projection_path is not None:
            p_to_data = p.projection_path
            num_s = p.file_name.find("{")
            num_e = p.file_name.find("}") + 1
            f_name = p.file_name.replace(p.file_name[num_s:num_e], "")
            proj_name = f_name.replace(p.im_ext, "{:s}Projection.tif")
            im_xy = imread(p_to_data + proj_name.format("xy")).transpose(
                2, 1, 0
            )
            im_xz = imread(p_to_data + proj_name.format("xz")).transpose(
                2, 1, 0
            )
            im_yz = imread(p_to_data + proj_name.format("yz")).transpose(
                2, 1, 0
            )
            self.viewer.add_image(
                im_xy,
                name=proj_name.format("xy"),
                scale=(1, p.voxel_size_out[1], p.voxel_size_out[0]),
            )
            self.viewer.add_image(
                im_yz,
                name=proj_name.format("xz"),
                scale=(1, p.voxel_size_out[2], p.voxel_size_out[1]),
                visible=False,
            )
            self.viewer.add_image(
                im_xz,
                name=proj_name.format("yz"),
                scale=(1, p.voxel_size_out[2], p.voxel_size_out[0]),
                visible=False,
            )

    def make_manual_parameterization(self):
        ### Paths definition
        path_data = self.make_file_search(
            label="Path to images", filter="", where="path_data_file", mode="d"
        )
        file_name = self.make_text_edit(
            label="File format",
            where="file_format",
            default="image{t:03d}.tiff",
        )
        trsf_folder = self.make_file_search(
            label="Trsf folder", filter="", where="trsf_folder", mode="d"
        )
        output_path = self.make_file_search(
            label="Output Format",
            filter="",
            where="output_path",
            default="_registered",
            mode="w",
        )
        projection_path = self.make_file_search(
            label="Path to projection",
            filter="",
            where="projection_path",
            mode="d",
        )
        path_tab = widgets.Container(
            widgets=[
                path_data,
                file_name,
                trsf_folder,
                output_path,
                projection_path,
            ],
            labels=False,
        )
        ### Geometry definition
        vox_res = self.make_voxel_label("Voxel size", "{axis}_vox")
        vox_res_out = self.make_voxel_label("Voxel size out", "{axis}_vox_out")
        geometry_tab = widgets.Container(
            widgets=[vox_res, vox_res_out], labels=False
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

        time_tab = widgets.Container(widgets=[f_tp, l_tp, r_tp], labels=False)

        ### Trsf parameterization
        trsf_type_label = widgets.Label(value="Transformation type")
        self.trsf_type = widgets.ComboBox(
            value="rigid", choices=["rigid", "affine"]
        )
        trsf_type = widgets.Container(
            widgets=[trsf_type_label, self.trsf_type], labels=False
        )

        compute_trsf = self.make_tick_box(
            "Compute Transformations", True, where="compute_trsf"
        )
        padding = self.make_tick_box("Padding", True, where="padding")
        recompute = self.make_tick_box(
            "Force recompute", True, where="recompute"
        )
        apply_trsf = self.make_tick_box(
            "Apply Trsf", True, where="apply_trsfs"
        )

        trsf_tab = widgets.Container(
            widgets=[compute_trsf, trsf_type, padding, apply_trsf, recompute],
            labels=False,
        )

        ### Advanced parameterization
        no_tp = self.make_text_edit(
            "Time points to skip:", "no_tp", "[]", eval_val=True
        )
        pre2D = self.make_tick_box("Pre 2D registration", False, where="pre2D")
        low_th = self.make_tick_box(
            "Low threshold", False, where="low_th_bool"
        )
        low_th_value_label = widgets.Label(value="Low threshold value:")
        self.low_th_val = widgets.IntText(value=250)
        low_th_value = widgets.Container(
            widgets=[low_th_value_label, self.low_th_val],
            layout="horizontal",
            labels=False,
        )
        reg_depth_label_start = widgets.Label(
            value="Registration depth start:"
        )
        self.reg_depth_start = widgets.IntText(value=6, min=2, max=6)
        r_depth_start = widgets.Container(
            widgets=[reg_depth_label_start, self.reg_depth_start],
            layout="horizontal",
            labels=False,
        )
        reg_depth_label_end = widgets.Label(value="Registration depth end:")
        self.reg_depth_end = widgets.IntText(value=3, min=1, max=5)
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
            labels=False,
        )

        plot_trsf = self.make_tick_box(
            "Plot trsf before apply", False, where="plot_trsf"
        )
        save_manual_setup = self.make_button(
            "Save parameters", self._save_json
        )

        clear_tick_manual = self.make_tick_box("Clear\nImages", default=False, where="clear_manual")
        start_reg_manual = widgets.Container(
            widgets=[
                clear_tick_manual,
                self.make_button("Run!", self._on_click_manual)
            ],
            labels=False,
            layout="horizontal"
        )
        
        sub_tab = QTabWidget()
        sub_tab.addTab(path_tab.native, "Paths")
        sub_tab.addTab(trsf_tab.native, "Trsf")
        sub_tab.addTab(geometry_tab.native, "Geometry")
        sub_tab.addTab(time_tab.native, "Time")
        sub_tab.addTab(advanced_tab.native, "Advanced")
        sub_tab.native = sub_tab
        sub_tab.name = ""
        manual_controler = widgets.Container(
            widgets=[sub_tab, plot_trsf, save_manual_setup, start_reg_manual],
            labels=False,
        )

        return manual_controler

    def _check_binaries(self):
        binary_list = [
            "applyTrsf",
            "blockmatching",
            "composeTrsf",
            "changeMultipleTrsfs",
        ]
        missed_binaries = []
        for bin in binary_list:
            if which(bin) is None:
                missed_binaries.append(bin)
        return missed_binaries

    def __init__(self, napari_viewer):
        super().__init__(napari_viewer)


class SpatialRegistrationWidget(RegistrationWidget):
    @property
    def link_to_parameters(self):
        self._link_to_parameters = {
            "path_to_data": self.path_data_file,
            "ref_im": self.ref_im,
            "flo_ims": self.flo_ims,
            "init_trsfs": self.init_trsfs,
            "trsf_types": self.trsf_types,
            "trsf_paths": self.trsf_paths,
            "ref_voxel": [self.x_vox, self.y_vox, self.z_vox],
            "flo_voxels": self.flo_voxels,
            "flo_im_sizes": None,
            "out_voxel": [self.x_out_vox, self.y_out_vox, self.z_out_vox],
            "out_pattern": self.out_pattern,
            "path_to_bin": "",
            "registration_depth_end": self.registration_depth_end,
            "registration_depth_start": self.registration_depth_start,
            "init_trsf_real_unit": True,
            "image_interpolation": self.image_interpolation,
            "apply_trsf": self.apply_trsf,
            "compute_trsf": self.compute_trsf,
            "test_init": self.test_init,
            "begin": 1,
            "end": 1,
            "time_tag": None,
            "bdv_unit": "microns",
            "do_bdv": False,
            "copy_ref": self.copy_ref,
            "bbox_out": self.bbox_out,
            "recompute": True,
        }
        return self._link_to_parameters

    @property
    def trsf_types(self):
        self._trsf_types = []
        if self.rigid.value:
            self._trsf_types.append("rigid")
        if self.affine.value:
            self._trsf_types.append("affine")
        if self.vectorfield.value:
            self._trsf_types.append("vectorfield")
        return self._trsf_types

    @trsf_types.setter
    def trsf_types(self, value):
        self._trsf_types = value
        self.rigid.value = False
        self.affine.value = False
        self.vectorfield.value = False
        for v in value:
            if isinstance(v, str):
                trsf = v.lower().strip()
                if trsf == "rigid":
                    self.rigid.value = True
                elif trsf == "affine":
                    self.affine.value = True
                elif trsf == "vectorfield":
                    self.vectorfield.value = True

    @property
    def flo_ims(self):
        self._flo_ims = []
        for i in range(5):
            if self.__dict__[f"to_use_{i+1}"].value:
                self._flo_ims.append(
                    self.__dict__[f"flo_im_{i+1}".lower()].value
                )
        return self._flo_ims

    @flo_ims.setter
    def flo_ims(self, value):
        self._flo_ims = value
        for angle_number, v in enumerate(value):
            self.__dict__[f"to_use_{angle_number+1}"].value = True
            self.__dict__[f"flo_im_{angle_number+1}".lower()].value = v

    @property
    def flo_voxels(self):
        self._flo_voxels = []
        for i in range(5):
            if self.__dict__[f"to_use_{i+1}"].value:
                flo_vox_tmp = []
                for axis in ["x", "y", "z"]:
                    flo_vox_tmp.append(
                        self.__dict__[f"{axis}_vox_flo_{i+1}"].value
                    )
                self._flo_voxels.append(flo_vox_tmp)
        return self._flo_voxels

    @flo_voxels.setter
    def flo_voxels(self, value):
        self._flo_voxels = value
        for angle_number, v in enumerate(value):
            self.__dict__[f"to_use_{angle_number+1}"].value = True
            for axis_pos, axis in enumerate(["x", "y", "z"]):
                self.__dict__[f"{axis}_vox_flo_{angle_number+1}"].value = v[
                    axis_pos
                ]

    @property
    def init_trsfs(self):
        self._init_trsfs = []
        for i in range(5):
            init_trsf = []
            if self.__dict__[f"to_use_{i+1}"].value:
                # Flip
                for axis in ["X", "Y", "Z"]:
                    if self.__dict__[f"{axis}_flip_{i+1}".lower()].value:
                        init_trsf.extend(["flip", f"{axis}"])
                # Rotation
                for axis in ["X", "Y", "Z"]:
                    if self.__dict__[f"rot_{axis}_{i+1}".lower()].value != 0:
                        init_trsf.extend(
                            [
                                "rot",
                                f"{axis}",
                                self.__dict__[
                                    f"rot_{axis}_{i+1}".lower()
                                ].value,
                            ]
                        )
                # Translation
                for axis in ["X", "Y", "Z"]:
                    if self.__dict__[f"trans_{axis}_{i+1}".lower()].value != 0:
                        init_trsf.extend(
                            [
                                "trans",
                                f"{axis}",
                                self.__dict__[
                                    f"trans_{axis}_{i+1}".lower()
                                ].value,
                            ]
                        )
                self._init_trsfs.append(init_trsf)
        return self._init_trsfs

    @init_trsfs.setter
    def init_trsfs(self, value):
        self._init_trsfs = value
        if not isinstance(value[0], list):
            value = [value]
        for axis in ["X", "Y", "Z"]:
            for i in range(5):
                self.__dict__[f"{axis}_flip_{i+1}".lower()].value = False
                self.__dict__[f"trans_{axis}_{i+1}".lower()].value = 0
                self.__dict__[f"rot_{axis}_{i+1}".lower()].value = 0
        for i, trsfs in enumerate(value):
            curr = 0
            while curr < len(trsfs):
                if trsfs[curr] == "flip":
                    self.__dict__[
                        f"{trsfs[curr+1].upper()}_flip_{i+1}".lower()
                    ].value = True
                    curr += 2
                else:
                    self.__dict__[
                        f"{trsfs[curr].lower()}_{trsfs[curr+1].upper()}_{i+1}".lower()
                    ].value = trsfs[curr + 2]
                    curr += 3

    def _on_click_manual(self):
        if self.clear_manual:
            for l in self.viewer.layers[:]:
                self.viewer.layers.remove(l)
        params = self.parameters
        nb_angles = 0
        for i in range(5):
            if self.__dict__[f"to_use_{i+1}"].value:
                nb_angles += 1
        params["trsf_paths"] = [
            params["trsf_paths"],
        ] * nb_angles
        tr = SpatialRegistration(params)
        print(self.init_trsfs)
        tr.run_trsf()
        p = tr.params[0]
        ref = imread(p.ref_out)
        vox = [v.value for v in self.link_to_parameters["out_voxel"]]
        color_maps = [
            "magenta",
            "cyan",
            "bop orange",
            "bop blue",
            "bop purple",
        ]
        self.viewer.add_image(
            ref,
            colormap="gray",
            scale=vox,
            name="Reference",
            blending="additive",
        )
        flos = []
        for i, p_flo in enumerate(p.flo_outs):
            im = imread(p_flo)
            flos.append(im)
            self.viewer.add_image(
                im,
                colormap=color_maps[i % len(color_maps)],
                scale=vox,
                name=f"Floating {i+1}",
                blending="additive",
            )

    def _on_click_file(self):
        if self.clear_file:
            for l in self.viewer.layers[:]:
                self.viewer.layers.remove(l)
        tr = SpatialRegistration(self.json_file.value)
        tr.run_trsf()
        p = tr.params[0]
        ref = imread(p.ref_out)
        vox = [v.value for v in self.link_to_parameters["out_voxel"]]
        color_maps = [
            "magenta",
            "cyan",
            "bop orange",
            "bop blue",
            "bop purple",
        ]
        self.viewer.add_image(
            ref,
            colormap="gray",
            scale=vox,
            name=f"Reference",
            blending="additive",
        )
        flos = []
        for i, p_flo in enumerate(p.flo_outs):
            im = imread(p_flo)
            flos.append(im)
            self.viewer.add_image(
                im,
                colormap=color_maps[i % len(color_maps)],
                scale=vox,
                name=f"Floating {i+1}",
                blending="additive",
            )

    def make_trans(self, label, where, min, max):
        trans_label = widgets.Label(value=label)
        trans_x_label = widgets.Label(value="X")
        self.__dict__[where.format(axis="x")] = widgets.FloatText(
            value=0, min=min, max=max
        )
        trans_x = widgets.Container(
            widgets=[trans_x_label, self.__dict__[where.format(axis="x")]],
            layout="horizontal",
            labels=False,
        )

        trans_y_label = widgets.Label(value="Y")
        self.__dict__[where.format(axis="y")] = widgets.FloatText(
            value=0, min=min, max=max
        )
        trans_y = widgets.Container(
            widgets=[trans_y_label, self.__dict__[where.format(axis="y")]],
            layout="horizontal",
            labels=False,
        )

        trans_z_label = widgets.Label(value="Z")
        self.__dict__[where.format(axis="z")] = widgets.FloatText(
            value=0, min=min, max=max
        )
        trans_z = widgets.Container(
            widgets=[trans_z_label, self.__dict__[where.format(axis="z")]],
            layout="horizontal",
            labels=False,
        )
        return widgets.Container(
            widgets=[trans_label, trans_x, trans_y, trans_z],
            layout="vertical",
            labels=False,
        )

    def make_angle_tab(self, angle_index):
        to_use = self.make_tick_box(
            "Use this angle",
            default=angle_index <= 1,
            where=f"to_use_{angle_index}",
        )
        to_use.native.layout().setSpacing(0)
        to_use.native.layout().setContentsMargins(0, 0, 0, 0)
        flo_vox = self.make_voxel_label(
            "Voxel size", "{axis}_vox_flo_" + f"{angle_index}"
        )
        flo_vox.native.layout().setSpacing(0)
        flo_vox.native.layout().setContentsMargins(0, 0, 0, 0)
        flo_im = self.make_text_edit(
            label=f"Floating image {angle_index}",
            where=f"flo_im_{angle_index}",
            default=f"image_{angle_index}.tiff",
        )
        flo_im.native.layout().setSpacing(0)
        flo_im.native.layout().setContentsMargins(0, 0, 0, 0)

        flip_label = widgets.Label(value="Flip:")
        x_flip = self.make_tick_box(
            "X", where=f"x_flip_{angle_index}", default=False
        )
        y_flip = self.make_tick_box(
            "Y", where=f"y_flip_{angle_index}", default=False
        )
        z_flip = self.make_tick_box(
            "Z", where=f"z_flip_{angle_index}", default=False
        )
        flip_tick = widgets.Container(
            widgets=[flip_label, x_flip, y_flip, z_flip],
            layout="horizontal",
            labels=False,
        )
        flip_tick.native.layout().setSpacing(0)
        flip_tick.native.layout().setContentsMargins(0, 0, 0, 0)
        # flip = widgets.Container(widgets=[flip_label, flip_tick], layout="vertical", labels=False)
        trans = self.make_trans(
            "Translation", "trans_{axis}_" + f"{angle_index}", -1e6, 1e6
        )
        trans.native.layout().setSpacing(0)
        trans.native.layout().setContentsMargins(0, 0, 0, 0)
        rot = self.make_trans(
            "Rotation (deg)", "rot_{axis}_" + f"{angle_index}", -360, 360
        )
        rot.native.layout().setSpacing(0)
        rot.native.layout().setContentsMargins(0, 0, 0, 0)
        basic = widgets.Container(
            widgets=[to_use, flo_im, flo_vox], labels=False
        )
        init_trsf = widgets.Container(widgets=[flip_tick, trans, rot])
        angle_i = QTabWidget()
        basic.native.layout().setSpacing(0)
        basic.native.layout().setContentsMargins(0, 0, 0, 0)
        init_trsf.native.layout().setSpacing(0)
        init_trsf.native.layout().setContentsMargins(0, 0, 0, 0)
        angle_i.addTab(basic.native, "Basic")
        angle_i.addTab(init_trsf.native, "Init trsf")
        angle_i.native = angle_i
        return angle_i

    def make_manual_parameterization(self):
        main_combobox = QComboBox()
        main_combobox.addItem("Global parameters")
        main_combobox.addItem("Reference")
        for i in range(5):
            main_combobox.addItem(f"Angle {i+1}")
        main_combobox._explicitly_hidden = False
        main_combobox.native = main_combobox

        main_stack = QStackedWidget()
        main_stack.native = main_stack

        # Global parameters
        ## Path
        path_data = self.make_file_search(
            label="Common path to images", filter=None, where="path_data_file"
        )

        trsf_paths = self.make_file_search(
            label="Path to save transformations",
            filter=None,
            where="trsf_paths",
            default="",
            mode="d",
        )

        out_pattern = self.make_text_edit(
            label="Output Pattern",
            where="out_pattern",
            default="_registered",
        )
        global_path_tab = widgets.Container(
            widgets=[path_data, trsf_paths, out_pattern]
        )

        ## Transformations types
        trsf_types_label = widgets.Label(value="Transformations types")
        rigid = self.make_tick_box("Rigid", True, where="rigid")
        affine = self.make_tick_box("Affine", True, where="affine")
        vectorfield = self.make_tick_box(
            "Non Linear", False, where="vectorfield"
        )
        out_voxel = self.make_voxel_label(
            "Output voxel size", "{axis}_out_vox"
        )
        global_trsf_tab = widgets.Container(
            widgets=[trsf_types_label, rigid, affine, vectorfield, out_voxel],
            labels=False,
        )

        ## Advanced parameters
        registration_depth_start_label = widgets.Label(
            value="Registration depth start:"
        )
        self.registration_depth_start = widgets.IntText(value=6, min=2, max=6)
        registration_depth_start = widgets.Container(
            widgets=[
                registration_depth_start_label,
                self.registration_depth_start,
            ],
            layout="horizontal",
            labels=False,
        )
        registration_depth_end_label = widgets.Label(
            value="Registration depth end:"
        )
        self.registration_depth_end = widgets.IntText(value=3, min=1, max=5)
        registration_depth_end = widgets.Container(
            widgets=[
                registration_depth_end_label,
                self.registration_depth_end,
            ],
            layout="horizontal",
            labels=False,
        )
        init_trsf_real_unit = self.make_tick_box(
            "Initial transformation in real units",
            True,
            where="init_trsf_real_unit",
        )

        image_interpolation_label = widgets.Label(value="Interpolation type")
        self.image_interpolation = widgets.ComboBox(
            value="linear", choices=["linear", "cspline"]
        )
        image_interpolation = widgets.Container(
            widgets=[image_interpolation_label, self.image_interpolation],
            labels=False,
        )

        apply_trsf = self.make_tick_box(
            "Apply transformation", True, "apply_trsf"
        )

        compute_trsf = self.make_tick_box(
            "Compute transformations", True, "compute_trsf"
        )

        copy_ref = self.make_tick_box(
            "Make a copy of the reference image", False, "copy_ref"
        )

        bbox_out = self.make_tick_box(
            "Build a fully englobing resulting frame", False, "bbox_out"
        )

        global_advanced_tab = widgets.Container(
            widgets=[
                registration_depth_start,
                registration_depth_end,
                init_trsf_real_unit,
                image_interpolation,
                apply_trsf,
                compute_trsf,
                copy_ref,
                bbox_out,
            ],
            labels=False,
        )
        global_advanced_tab.native.layout().setSpacing(0)
        global_advanced_tab.native.layout().setContentsMargins(0, 0, 0, 0)

        global_path_tab.native.layout().setSpacing(0)
        global_path_tab.native.layout().setContentsMargins(0, 0, 0, 0)

        global_trsf_tab.native.layout().setSpacing(0)
        global_trsf_tab.native.layout().setContentsMargins(0, 0, 0, 0)

        global_tab = QTabWidget()
        global_tab.addTab(global_path_tab.native, "Path")
        global_tab.addTab(global_trsf_tab.native, "Trsf types")
        global_tab.addTab(global_advanced_tab.native, "Advanced")

        # Reference
        ref_im = self.make_text_edit(
            label="Reference image name",
            where="ref_im",
            default="image_0.tiff",
        )

        ref_voxel = self.make_voxel_label(
            "Reference voxel size", where="{axis}_vox"
        )

        ref_im.native.layout().setSpacing(0)
        ref_im.native.layout().setContentsMargins(0, 0, 0, 0)
        ref_voxel.native.layout().setSpacing(0)
        ref_voxel.native.layout().setContentsMargins(0, 0, 0, 0)

        ref_tab = widgets.Container(widgets=[ref_im, ref_voxel])
        ref_tab.native.layout().setSpacing(0)
        ref_tab.native.layout().setContentsMargins(0, 0, 0, 0)

        main_stack.addWidget(global_tab)
        main_stack.addWidget(ref_tab.native)

        # Floating images
        for i in range(5):
            angle = self.make_angle_tab(i + 1)
            main_stack.addWidget(angle.native)

        save_manual_setup = self.make_button(
            "Save parameters", self._save_json
        )
        test_init = self.make_tick_box(
            "Only test initial transformation", False, "test_init"
        )

        clear_tick_manual = self.make_tick_box("Clear\nImages", default=False, where="clear_manual")
        start_reg_manual = widgets.Container(
            widgets=[
                clear_tick_manual,
                self.make_button("Run!", self._on_click_manual)
            ],
            labels=False,
            layout="horizontal"
        )
        
        # start_reg_manual = self.make_button("Run!", self._on_click_manual)

        main_combobox.currentIndexChanged.connect(main_stack.setCurrentIndex)
        main_combobox.name = "main_combobox"
        main_stack.name = "main_stack"
        main_control = widgets.Container(
            widgets=[
                main_combobox,
                main_stack,
            ],
            labels=False,
        )

        manual_controler = widgets.Container(
            widgets=[
                main_control,
                test_init,
                save_manual_setup,
                start_reg_manual,
            ],
            labels=False,
        )

        return manual_controler

    def _check_binaries(self):
        binary_list = [
            "applyTrsf",
            "blockmatching",
            "changeMultipleTrsfs",
            "invTrsf",
        ]
        missed_binaries = []
        for bin in binary_list:
            if which(bin) is None:
                missed_binaries.append(bin)
        return missed_binaries

    def __init__(self, napari_viewer):
        super().__init__(napari_viewer)
