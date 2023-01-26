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
)
from registrationtools import SpatialRegistration, TimeRegistration
from pathlib import Path
from IO import imread
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
            if k in ["voxel_size", "voxel_size_out"]:
                val = []
                for vs in v:
                    val.append(vs.value)
            else:
                if hasattr(v, 'value'):
                    val = v.value
                else:
                    val = v

                if isinstance(val, Path):
                    val = str(val)
            self._parameters[k] = val
        return self._parameters

    def make_file_search(self, label, filter, where=None, default=None, mode='r'):
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

    def make_voxel_label(self, label="Voxel size", where='{axis}_vox'):
        vox_label = widgets.Label(value=label)
        x_label = widgets.Label(value="x")
        self.__dict__[where.format(axis="x")] = widgets.FloatText(value=1)
        cont_x = widgets.Container(
            widgets=[
                x_label,
                self.__dict__[where.format(axis="x")]
            ],
            layout="horizontal",
            labels=False,
        )
        y_label = widgets.Label(value="y")
        self.__dict__[where.format(axis="y")] = widgets.FloatText(value=1)
        cont_y = widgets.Container(
            widgets=[
                y_label,
                self.__dict__[where.format(axis="y")]
            ],
            layout="horizontal",
            labels=False,
        )
        z_label = widgets.Label(value="z")
        self.__dict__[where.format(axis="z")] = widgets.FloatText(value=6)
        cont_z = widgets.Container(
            widgets=[
                z_label,
                self.__dict__[where.format(axis="z")]
            ],
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
        for k, v in json_params.items():
            if k in self.link_to_parameters and hasattr(
                self.link_to_parameters[k], "value"
            ):
                self.link_to_parameters[k].value = v
            elif k in ["voxel_size", "voxel_size_out"]:
                for i, vs in enumerate(v):
                    self.link_to_parameters[k][i].value = vs
            elif k == "trsf_types":
                self.trsf_types = v
                self.link_to_parameters[k] = self.trsf_types

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
                json.dump(self.parameters, f)

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

    def __init__(self, napari_viewer):
        super().__init__()
        self.viewer = napari_viewer

        ## json file parameterization:
        json = self.make_file_search("Parameter file", "*.json", "json_file")
        json.changed.connect(self._propagate_json)
        start_reg = self.make_button("Run!", self._on_click_file)
        self.file_params = widgets.Container(
            widgets=[json, start_reg], labels=False
        )
        self.file_params.native.layout().addStretch(1)

        ## Manual parameterization
        manual_controler = self.make_manual_parameterization()
        
        self.tab_controls = QTabWidget()
        self.tab_controls.addTab(self.file_params.native, "From File")
        self.tab_controls.addTab(manual_controler.native, "Manual")

        layout = QVBoxLayout()
        layout.addStretch(1)
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
        }
        return self._link_to_parameters

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

    def make_manual_parameterization(self):
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
                projection_path,
            ],
            labels=False,
        )
        ### Geometry definition
        vox_res = self.make_voxel_label("Voxel size", "{axis}_vox")
        vox_res_out = self.make_voxel_label("Voxel size out", "{axis}_vox_out")
        geometry_tab = widgets.Container(widgets=[vox_res, vox_res_out], labels=False)

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
        no_tp = self.make_text_edit("Time points to skip:", "no_tp", "[]", eval_val=True)
        pre2D = self.make_tick_box("Pre 2D registration", True, where="pre2D")
        low_th = self.make_tick_box("Low threshold", True, where="low_th")
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
            "Plot trsf before apply", True, where="plot_trsf"
        )
        save_manual_setup = self.make_button(
            "Save parameters", self._save_json
        )
        start_reg_manual = self.make_button("Run!", self._on_click_manual)

        sub_tab = QTabWidget()
        sub_tab.addTab(path_tab.native, "Paths")
        sub_tab.addTab(trsf_tab.native, "Trsf")
        sub_tab.addTab(geometry_tab.native, "Geometry")
        sub_tab.addTab(time_tab.native, "Time")
        sub_tab.addTab(advanced_tab.native, "Advanced")
        sub_tab.native = sub_tab
        manual_controler = widgets.Container(
            widgets=[sub_tab, plot_trsf, save_manual_setup, start_reg_manual],
            labels=False,
        )

        return manual_controler

    def __init__(self, napari_viewer):
        super().__init__(napari_viewer)

class SpatialRegistrationWidget(RegistrationWidget):
    def split_string(self, string):
        sep = self.sep.value
        return [k.strip() for k in string.split(sep)]

    @property
    def link_to_parameters(self):
        self._link_to_parameters = {
            "path_to_data": self.path_data_file,
            "ref_im": self.ref_im,
            "flo_ims": self.flo_ims,
            "init_trsfs": self.init_trsfs,
            "trsf_types": self.trsf_types,
            "trsf_paths": self.trsf_paths,
            "ref_voxel": self.ref_voxel,
            "flo_voxels": self.flo_voxels,
            "flo_im_sizes": None,
            "out_voxel": self.out_voxel,
            "out_pattern": self.out_pattern,
            "path_to_bin": "",
            "registration_depth": self.registration_depth,
            "init_trsf_real_unit": True,
            "image_interpolation": self.image_interpolation,
            "apply_trsf": self.apply_trsf,
            "compute_trsf": self.compute_trsf,
            "test_init": self.test_init,
            "begin": 1,
            "end": 1,
            "time_tag": "TM",
            "bdv_unit": "microns",
            "do_bdv": False,
            "copy_ref": self.copy_ref,
            "bbox_out": self.bbox_out
        }
        return self._link_to_parameters

    @property
    def trsf_types(self):
        self._trsf_types = []
        if self.rigid.value:
            self._trsf_types.append('rigid')
        if self.affine.value:
            self._trsf_types.append('affine')
        if self.vectorfield.value:
            self._trsf_types.append('vectorfield')
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
                if trsf == 'rigid':
                    self.rigid.value = True
                elif trsf == 'affine':
                    self.affine.value = True
                elif trsf == 'vectorfield':
                    self.vectorfield.value = True

    def _on_click_manual(self):
        params = self.parameters
        params["trsf_paths"] = [params["trsf_paths"]]
        tr = SpatialRegistration(params)
        tr.run_trsf()
        p = tr.params[0]
        ref = imread(p.ref_out)
        vox = self.out_voxel.value
        color_maps = ['magenta', 'cyan', 'bop blue', 'bop orange', 'bop purple']
        self.viewer.add_image(
            ref,
            colormap='gray',
            scale=vox,
            name="Reference",
            opacity=.5
        )
        flos = []
        for i, p_flo in enumerate(p.flo_outs):
            im = imread(p_flo)
            flos.append(im)
            self.viewer.add_image(
                im,
                colormap=color_maps[i%len(color_maps)],
                scale=vox,
                name=f"Floating {i}",
                opacity=.5
            )

    def _on_click_file(self):
        tr = SpatialRegistration(self.json_file.value)
        tr.run_trsf()
        p = tr.params[0]
        ref = imread(p.ref_out)
        vox = self.out_voxel.value
        color_maps = ['magenta', 'cyan', 'bop blue', 'bop orange', 'bop purple']
        self.viewer.add_image(
            ref,
            colormap='gray',
            scale=vox,
            name=f"Reference",
            opacity=.5
        )
        flos = []
        for i, p_flo in enumerate(p.flo_outs):
            im = imread(p_flo)
            flos.append(im)
            self.viewer.add_image(
                im,
                colormap=color_maps[i%len(color_maps)],
                scale=vox,
                name=f"Flo {i}",
                opacity=.5
            )


    def make_manual_parameterization(self):

        ### Paths definition
        path_data = self.make_file_search(
            label="Path to images", filter=None, where="path_data_file"
        )
        ref_im = self.make_text_edit(
            label="Reference angle",
            where="ref_im",
            default="image_0.tiff",
        )
        flo_ims = self.make_text_edit(
            label="Floating images",
            where="flo_ims",
            default='["image_1.tiff", "image_2.tiff", "image_3.tiff"]',
            eval_val=True
        )
        trsf_paths = self.make_file_search(
            label="Path to save transformations",
            filter=None,
            where="trsf_paths",
            default="",
            mode='d'
        )
        
        out_pattern = self.make_text_edit(
            label="Output Pattern",
            where="out_pattern",
            default="_registered",
        )
        path_tab = widgets.Container(
            widgets=[
                path_data,
                ref_im,
                flo_ims,
                trsf_paths,
                out_pattern,
            ],
            labels=False,
        )

        ### Transformations
        init_trsfs = self.make_text_edit(
            label="Initial transformations",
            where="init_trsfs",
            default=('[["flip", "Y", "flip", "Z", "trans", "Z", -72],\n'
                     ' ["flip", "Y"],\n'
                     ' ["flip", "Z"]]'),
            eval_val=True
        )

        trsf_types_label = widgets.Label(value="Transformations types")
        rigid = self.make_tick_box(
            "Rigid",
            True,
            where="rigid"
        )
        affine = self.make_tick_box(
            "Affine",
            True,
            where="affine"
        )
        vectorfield = self.make_tick_box(
            "Non Linear",
            False,
            where="vectorfield"
        )
        trsf_types = widgets.Container(
            widgets=[
                trsf_types_label,
                rigid,
                affine,
                vectorfield
            ]
        )

        trsf_folder = self.make_file_search(
            label="Trsf folder", filter=None, where="trsf_folder", mode="d"
        )
        trsf_tab = widgets.Container(
            widgets=[
                init_trsfs,
                trsf_types,
                trsf_folder
            ],
            labels=False,
        )

        ### Geometry definition
        ref_voxel = self.make_text_edit(
            label="Reference voxel",
            where="ref_voxel",
            default='.5, .5, 1',
            eval_val=True
        )

        flo_voxels = self.make_text_edit(
            label="Floating voxel",
            where="flo_voxels",
            default='[.5, .5, 1], [.5, .5, 1], [.2, .2, 1]',
            eval_val=True
        )

        out_voxel = self.make_text_edit(
            label="Output voxel",
            where="out_voxel",
            default='1, 1, 1',
            eval_val=True
        )

        geometry_tab = widgets.Container(
            widgets=[
                ref_voxel,
                flo_voxels,
                out_voxel
            ],
            labels=False,
        )

        ### Advancedlink_to_parameters
        registration_depth_label = widgets.Label(
            value="Registration depth start:"
        )
        self.registration_depth = widgets.IntText(value=3, min=1, max=5)
        registration_depth = widgets.Container(
            widgets=[registration_depth_label, self.registration_depth],
            layout="horizontal",
            labels=False,
        )
        init_trsf_real_unit = self.make_tick_box(
            "Initial transformation in real units",
            True,
            where="init_trsf_real_unit"
        )

        image_interpolation_label = widgets.Label(value="Interpolation type")
        self.image_interpolation = widgets.ComboBox(
            value="linear", choices=["linear", "cspline"]
        )
        image_interpolation = widgets.Container(
            widgets=[image_interpolation_label, self.image_interpolation], labels=False
        )

        apply_trsf = self.make_tick_box(
            "Apply transformation",
            True,
            "apply_trsf"
        )

        compute_trsf = self.make_tick_box(
            "Compute transformations",
            True,
            "compute_trsf"
        )

        copy_ref = self.make_tick_box(
            "Make a copy of the reference image",
            False,
            "copy_ref"
        )

        bbox_out = self.make_tick_box(
            "Build a fully englobing resulting frame",
            False,
            "bbox_out"
        )

        advanced_tab = widgets.Container(
            widgets=[
                registration_depth,
                init_trsf_real_unit,
                image_interpolation,
                apply_trsf,
                compute_trsf,
                copy_ref,
                bbox_out,
            ],
            labels=False,
        )

        save_manual_setup = self.make_button(
            "Save parameters", self._save_json
        )
        test_init = self.make_tick_box(
            "Only test initial transformation",
            False,
            "test_init"
        )
        start_reg_manual = self.make_button("Run!", self._on_click_manual)

        sub_tab = QTabWidget()
        sub_tab.addTab(path_tab.native, "Paths")
        sub_tab.addTab(trsf_tab.native, "Trsf")
        sub_tab.addTab(geometry_tab.native, "Geometry")
        sub_tab.addTab(advanced_tab.native, "Advanced")
        sub_tab.native = sub_tab
        manual_controler = widgets.Container(
            widgets=[sub_tab, test_init, save_manual_setup, start_reg_manual],
            labels=False,
        )

        return manual_controler

    def __init__(self, napari_viewer):
        super().__init__(napari_viewer)
