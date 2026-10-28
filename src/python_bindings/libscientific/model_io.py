"""
model_io.py - Python FFI bindings for LibScientific model container serialization.

Provides high-level Python interfaces for saving and importing standalone models (.lsm)
and multi-model project bundles (.lsb) with zero external dependencies.

Copyright (C) 2026 designed, written and maintained by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
"""

import ctypes
from enum import IntEnum
from typing import Any

from libscientific.loadlibrary import load_libscientific_library

lsci = load_libscientific_library()


class ModelType(IntEnum):
    """
    Enumeration of model types supported in the LSCM container format.
    """

    UNKNOWN = 0
    PCA = 1
    PLS = 2
    CPCA = 3
    MLR = 4
    LDA = 5
    ICA = 6
    UPCA = 7
    UPLS = 8
    ANY = 255


# LSCM status codes
LSCM_OK = 0
LSCM_ERR_FILE_NOT_FOUND = -1
LSCM_ERR_IO = -2
LSCM_ERR_INVALID_MAGIC = -3
LSCM_ERR_UNSUPPORTED_VER = -4
LSCM_ERR_INVALID_TYPE = -5
LSCM_ERR_TYPE_MISMATCH = -6
LSCM_ERR_NOT_FOUND = -7
LSCM_ERR_MEM = -8
LSCM_ERR_PARAM = -9
LSCM_ERR_CORRUPT = -10

_STATUS_MESSAGES = {
    0: "Success",
    -1: "File not found",
    -2: "I/O error during read/write",
    -3: "Invalid container magic (not an LSCM file)",
    -4: "Unsupported container format version",
    -5: "Invalid or unknown model type",
    -6: "Model type mismatch",
    -7: "Model identifier not found in bundle",
    -8: "Memory allocation failure",
    -9: "Invalid parameter or NULL pointer",
    -10: "Container TOC or payload corrupted",
}


def _check_status(status: int, action: str = "Operation"):
    if status != LSCM_OK:
        msg = _STATUS_MESSAGES.get(status, f"Unknown error code {status}")
        if status == LSCM_ERR_FILE_NOT_FOUND:
            raise FileNotFoundError(f"{action} failed: {msg}")
        elif status in (LSCM_ERR_PARAM, LSCM_ERR_INVALID_TYPE, LSCM_ERR_TYPE_MISMATCH):
            raise ValueError(f"{action} failed: {msg}")
        elif status == LSCM_ERR_MEM:
            raise MemoryError(f"{action} failed: {msg}")
        else:
            raise OSError(f"{action} failed: {msg}")


class MODEL_FILE_INFO(ctypes.Structure):
    """
    Metadata descriptor for inspecting files without full payload loading.
    """

    _fields_ = [
        ("num_models", ctypes.c_uint32),
        ("is_single", ctypes.c_bool),
        ("primary_type", ctypes.c_uint32),
        ("primary_name", ctypes.c_char * 64),
        ("format_version", ctypes.c_uint32),
        ("total_size", ctypes.c_uint64),
    ]


class MODEL_ENTRY(ctypes.Structure):
    """
    Single entry within a multi-model bundle.
    """

    _fields_ = [
        ("name", ctypes.c_char * 64),
        ("type", ctypes.c_uint32),
        ("model", ctypes.c_void_p),
    ]


class MODEL_BUNDLE(ctypes.Structure):
    """
    Multi-model bundle container.
    """

    _fields_ = [
        ("count", ctypes.c_size_t),
        ("capacity", ctypes.c_size_t),
        ("entries", ctypes.POINTER(MODEL_ENTRY)),
    ]


# C API Declarations
lsci.SaveModel.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_uint32, ctypes.c_void_p]
lsci.SaveModel.restype = ctypes.c_int

lsci.ImportModel.argtypes = [ctypes.c_char_p, ctypes.c_uint32, ctypes.c_void_p]
lsci.ImportModel.restype = ctypes.c_int

lsci.NewModelBundle.argtypes = []
lsci.NewModelBundle.restype = ctypes.POINTER(MODEL_BUNDLE)

lsci.ModelBundleAddModel.argtypes = [
    ctypes.POINTER(MODEL_BUNDLE),
    ctypes.c_char_p,
    ctypes.c_uint32,
    ctypes.c_void_p,
]
lsci.ModelBundleAddModel.restype = ctypes.c_int

lsci.ModelBundleGetModel.argtypes = [
    ctypes.POINTER(MODEL_BUNDLE),
    ctypes.c_char_p,
    ctypes.POINTER(ctypes.c_uint32),
]
lsci.ModelBundleGetModel.restype = ctypes.c_void_p

lsci.ModelBundleGetModelByIndex.argtypes = [
    ctypes.POINTER(MODEL_BUNDLE),
    ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_char_p),
    ctypes.POINTER(ctypes.c_uint32),
]
lsci.ModelBundleGetModelByIndex.restype = ctypes.c_void_p

lsci.DelModelBundle.argtypes = [ctypes.POINTER(ctypes.POINTER(MODEL_BUNDLE)), ctypes.c_int]
lsci.DelModelBundle.restype = None

lsci.SaveModelBundle.argtypes = [ctypes.c_char_p, ctypes.POINTER(MODEL_BUNDLE)]
lsci.SaveModelBundle.restype = ctypes.c_int

lsci.ImportModelBundle.argtypes = [ctypes.c_char_p]
lsci.ImportModelBundle.restype = ctypes.POINTER(MODEL_BUNDLE)

lsci.InspectModelFile.argtypes = [ctypes.c_char_p, ctypes.POINTER(MODEL_FILE_INFO)]
lsci.InspectModelFile.restype = ctypes.c_int

lsci.ImportModelAuto.argtypes = [
    ctypes.c_char_p,
    ctypes.POINTER(ctypes.c_void_p),
    ctypes.POINTER(ctypes.c_uint32),
    ctypes.POINTER(ctypes.POINTER(MODEL_BUNDLE)),
]
lsci.ImportModelAuto.restype = ctypes.c_int

lsci.SavePCAModel.argtypes = [ctypes.c_char_p, ctypes.c_void_p]
lsci.SavePCAModel.restype = ctypes.c_int

lsci.ImportPCAModel.argtypes = [ctypes.c_char_p, ctypes.c_void_p]
lsci.ImportPCAModel.restype = ctypes.c_int

lsci.SavePLSModel.argtypes = [ctypes.c_char_p, ctypes.c_void_p]
lsci.SavePLSModel.restype = ctypes.c_int

lsci.ImportPLSModel.argtypes = [ctypes.c_char_p, ctypes.c_void_p]
lsci.ImportPLSModel.restype = ctypes.c_int

lsci.SaveCPCAModel.argtypes = [ctypes.c_char_p, ctypes.c_void_p]
lsci.SaveCPCAModel.restype = ctypes.c_int

lsci.ImportCPCAModel.argtypes = [ctypes.c_char_p, ctypes.c_void_p]
lsci.ImportCPCAModel.restype = ctypes.c_int

lsci.ModelTypeToString.argtypes = [ctypes.c_uint32]
lsci.ModelTypeToString.restype = ctypes.c_char_p

lsci.ModelTypeFromString.argtypes = [ctypes.c_char_p]
lsci.ModelTypeFromString.restype = ctypes.c_uint32


# Low-level C function aliases matching C API naming convention
SaveModel = lsci.SaveModel
ImportModel = lsci.ImportModel
SaveModelBundle = lsci.SaveModelBundle
ImportModelBundle = lsci.ImportModelBundle
InspectModelFile = lsci.InspectModelFile
ImportModelAuto = lsci.ImportModelAuto
NewModelBundle = lsci.NewModelBundle
ModelBundleAddModel = lsci.ModelBundleAddModel
ModelBundleGetModel = lsci.ModelBundleGetModel
ModelBundleGetModelByIndex = lsci.ModelBundleGetModelByIndex
DelModelBundle = lsci.DelModelBundle
SavePCAModel = lsci.SavePCAModel
ImportPCAModel = lsci.ImportPCAModel
SavePLSModel = lsci.SavePLSModel
ImportPLSModel = lsci.ImportPLSModel
SaveCPCAModel = lsci.SaveCPCAModel
ImportCPCAModel = lsci.ImportCPCAModel


def _extract_c_model_and_type(
    obj: Any, explicit_type: ModelType | None = None
) -> tuple[ctypes.c_void_p, ModelType]:
    """
    Extracts the underlying C pointer and ModelType from a high-level wrapper or ctypes struct.
    """
    cls_name = type(obj).__name__
    if hasattr(obj, "model"):
        # High-level wrapper class (PCA, PLS, CPCA, etc.)
        if cls_name == "PCA":
            return ctypes.cast(obj.model, ctypes.c_void_p), ModelType.PCA
        elif cls_name == "PLS":
            return ctypes.cast(obj.model, ctypes.c_void_p), ModelType.PLS
        elif cls_name == "CPCA":
            return ctypes.cast(obj.model, ctypes.c_void_p), ModelType.CPCA

    if hasattr(obj, "_type_"):
        # Direct ctypes pointer
        inner_name = obj._type_.__name__
        if inner_name == "PCAMODEL":
            return ctypes.cast(obj, ctypes.c_void_p), ModelType.PCA
        elif inner_name == "PLSMODEL":
            return ctypes.cast(obj, ctypes.c_void_p), ModelType.PLS
        elif inner_name == "CPCAMODEL":
            return ctypes.cast(obj, ctypes.c_void_p), ModelType.CPCA

    if explicit_type is not None and explicit_type != ModelType.UNKNOWN:
        return ctypes.cast(obj, ctypes.c_void_p), explicit_type

    raise ValueError(f"Unable to automatically determine model type for {obj} ({cls_name})")


def _wrap_c_model(ptr_val: Any, model_type: ModelType) -> Any:
    """
    Wraps an allocated C model pointer into its corresponding Python class.
    """
    if not ptr_val:
        return None

    if model_type == ModelType.PCA:
        from libscientific.pca import PCA, PCAMODEL

        pca_ptr = ctypes.cast(ptr_val, ctypes.POINTER(PCAMODEL))
        pca_inst = PCA.__new__(PCA)
        pca_inst.model = pca_ptr
        pca_inst.scaling = 1
        if pca_ptr and bool(pca_ptr.contents.scores):
            pca_inst.npc = pca_ptr.contents.scores[0].col
        else:
            pca_inst.npc = 0
        return pca_inst
    elif model_type == ModelType.PLS:
        from libscientific.pls import PLS, PLSMODEL

        pls_ptr = ctypes.cast(ptr_val, ctypes.POINTER(PLSMODEL))
        pls_inst = PLS.__new__(PLS)
        pls_inst.model = pls_ptr
        pls_inst.xscaling = 1
        pls_inst.yscaling = 0
        if pls_ptr and bool(pls_ptr.contents.xscores):
            pls_inst.nlv = pls_ptr.contents.xscores[0].col
        else:
            pls_inst.nlv = 0
        return pls_inst
    elif model_type == ModelType.CPCA:
        from libscientific.cpca import CPCA, CPCAMODEL

        cpca_ptr = ctypes.cast(ptr_val, ctypes.POINTER(CPCAMODEL))
        cpca_inst = CPCA.__new__(CPCA)
        cpca_inst.model = cpca_ptr
        return cpca_inst
    else:
        return ptr_val


class ModelBundle:
    """
    Container representing a multi-model LibScientific project bundle (.lsb).
    Provides dictionary-like access to contained models.
    """

    def __init__(self, models: dict[str, Any] | list[Any] | tuple[Any, ...] | None = None):
        self._models: dict[str, Any] = {}
        if models:
            if isinstance(models, dict):
                for k, v in models.items():
                    self.add_model(k, v)
            elif isinstance(models, (list, tuple)):
                for idx, v in enumerate(models):
                    self.add_model(f"model_{idx}", v)

    def add_model(self, name: str, model: Any, model_type: ModelType | None = None):
        """
        Add a model to the bundle with a unique name identifier.
        """
        if not name or not isinstance(name, str):
            raise ValueError("Model name must be a non-empty string")
        self._models[name] = model

    def get_model(self, name: str) -> Any:
        """
        Retrieve a model by name.
        """
        return self._models[name]

    def __getitem__(self, name: str) -> Any:
        return self._models[name]

    def __setitem__(self, name: str, model: Any):
        self.add_model(name, model)

    def __delitem__(self, name: str):
        del self._models[name]

    def __contains__(self, name: str) -> bool:
        return name in self._models

    def __len__(self) -> int:
        return len(self._models)

    def __iter__(self):
        return iter(self._models)

    def keys(self):
        return self._models.keys()

    def values(self):
        return self._models.values()

    def items(self):
        return self._models.items()

    def save(self, filepath: str):
        """
        Save this bundle to an .lsb file.
        """
        save_model_bundle(filepath, self)

    def __repr__(self) -> str:
        names = ", ".join(repr(k) for k in self._models.keys())
        return f"<ModelBundle count={len(self._models)} models=[{names}]>"


def inspect_model_file(filepath: str) -> dict[str, Any]:
    """
    Inspect a model file (.lsm or .lsb) header without loading matrix/tensor payloads.

    Returns:
        dict containing:
            - num_models (int)
            - is_single (bool)
            - primary_type (ModelType)
            - primary_name (str)
            - format_version (int)
            - total_size (int)
    """
    info = MODEL_FILE_INFO()
    st = lsci.InspectModelFile(filepath.encode("utf-8"), ctypes.byref(info))
    _check_status(st, f"InspectModelFile('{filepath}')")
    return {
        "num_models": int(info.num_models),
        "is_single": bool(info.is_single),
        "primary_type": ModelType(info.primary_type),
        "primary_name": info.primary_name.decode("utf-8", errors="replace"),
        "format_version": int(info.format_version),
        "total_size": int(info.total_size),
    }


def is_bundle(filepath: str) -> bool:
    """
    Check if a model file is a multi-model bundle.
    """
    info = inspect_model_file(filepath)
    return not info["is_single"]


def save_model(
    filepath: str, model: Any, model_name: str = "default", model_type: ModelType | None = None
):
    """
    Save a single model into a standalone .lsm file.

    Parameters:
        filepath (str): Destination file path (e.g. 'model.lsm')
        model: Model object (e.g. PCA, PLS, CPCA instance or ctypes pointer)
        model_name (str): Identifier tag for the model
        model_type (ModelType, optional): Explicit model type if passing raw pointer
    """
    c_ptr, mtype = _extract_c_model_and_type(model, model_type)
    fp = filepath.encode("utf-8")
    nm = model_name.encode("utf-8")
    st = lsci.SaveModel(fp, nm, int(mtype), c_ptr)
    _check_status(st, f"SaveModel('{filepath}')")


def import_model(
    filepath: str, expected_type: ModelType = ModelType.ANY, target_model: Any = None
) -> Any:
    """
    Import a single model from a file (.lsm or .lsb with 1 model).

    If target_model is supplied, its internal weights are overwritten with the loaded model.
    Otherwise, a new instance of the appropriate class (e.g. PCA, PLS) is automatically created.

    Parameters:
        filepath (str): Source file path
        expected_type (ModelType): Optional type constraint (default: ModelType.ANY)
        target_model: Optional existing model instance to populate in-place

    Returns:
        The populated model instance.
    """
    fp = filepath.encode("utf-8")
    if target_model is not None:
        c_ptr, mtype = _extract_c_model_and_type(target_model, expected_type)
        st = lsci.ImportModel(fp, int(mtype), c_ptr)
        _check_status(st, f"ImportModel('{filepath}')")
        if (
            hasattr(target_model, "npc")
            and hasattr(target_model, "model")
            and target_model.model.contents.scores
        ):
            target_model.npc = target_model.model.contents.scores[0].col
        elif (
            hasattr(target_model, "nlv")
            and hasattr(target_model, "model")
            and target_model.model.contents.xscores
        ):
            target_model.nlv = target_model.model.contents.xscores[0].col
        return target_model
    else:
        info = inspect_model_file(filepath)
        mtype = info["primary_type"]
        if expected_type != ModelType.ANY and mtype != expected_type:
            raise ValueError(
                f"Type mismatch: file contains {mtype.name}, expected {expected_type.name}"
            )

        if mtype == ModelType.PCA:
            from libscientific.pca import PCA

            pca_inst = PCA()
            st = lsci.ImportPCAModel(fp, pca_inst.model)
            _check_status(st, f"ImportPCAModel('{filepath}')")
            if pca_inst.model and pca_inst.model.contents.scores:
                pca_inst.npc = pca_inst.model.contents.scores[0].col
            return pca_inst
        elif mtype == ModelType.PLS:
            from libscientific.pls import PLS

            pls_inst = PLS()
            st = lsci.ImportPLSModel(fp, pls_inst.model)
            _check_status(st, f"ImportPLSModel('{filepath}')")
            if pls_inst.model and pls_inst.model.contents.xscores:
                pls_inst.nlv = pls_inst.model.contents.xscores[0].col
            return pls_inst
        elif mtype == ModelType.CPCA:
            from libscientific.cpca import CPCA

            cpca_inst = CPCA()
            st = lsci.ImportCPCAModel(fp, cpca_inst.model)
            _check_status(st, f"ImportCPCAModel('{filepath}')")
            return cpca_inst
        else:
            out_model = ctypes.c_void_p(None)
            out_type = ctypes.c_uint32(0)
            out_bundle = ctypes.POINTER(MODEL_BUNDLE)()
            st = lsci.ImportModelAuto(
                fp, ctypes.byref(out_model), ctypes.byref(out_type), ctypes.byref(out_bundle)
            )
            _check_status(st, f"ImportModelAuto('{filepath}')")
            return _wrap_c_model(out_model, ModelType(out_type.value))


def save_model_bundle(
    filepath: str, bundle: ModelBundle | dict[str, Any] | list[Any] | tuple[Any, ...]
):
    """
    Save multiple models into a single multi-model project file (.lsb).

    Parameters:
        filepath (str): Destination file path (e.g. 'project.lsb')
        bundle: ModelBundle, dictionary, or sequence of model instances
    """
    fp = filepath.encode("utf-8")
    c_bundle = lsci.NewModelBundle()
    if not c_bundle:
        raise MemoryError("Failed to allocate C MODEL_BUNDLE")

    try:
        if isinstance(bundle, ModelBundle):
            items = bundle.items()
        elif isinstance(bundle, dict):
            items = bundle.items()
        elif isinstance(bundle, (list, tuple)):
            items = [(f"model_{i}", m) for i, m in enumerate(bundle)]
        else:
            raise TypeError("Bundle must be a ModelBundle, dict, or list of models")

        for name, model in items:
            c_ptr, mtype = _extract_c_model_and_type(model)
            st = lsci.ModelBundleAddModel(c_bundle, name.encode("utf-8"), int(mtype), c_ptr)
            _check_status(st, f"ModelBundleAddModel('{name}')")

        st = lsci.SaveModelBundle(fp, c_bundle)
        _check_status(st, f"SaveModelBundle('{filepath}')")
    finally:
        # DelModelBundle with free_contained_models=0: leaves Python models alive
        lsci.DelModelBundle(ctypes.byref(c_bundle), 0)


def import_model_bundle(filepath: str) -> ModelBundle:
    """
    Import a multi-model project bundle from an .lsb file.

    Returns:
        ModelBundle mapping model names to model instances.
    """
    fp = filepath.encode("utf-8")
    c_bundle = lsci.ImportModelBundle(fp)
    if not c_bundle:
        raise OSError(f"ImportModelBundle failed to read '{filepath}'")

    bundle = ModelBundle()
    count = c_bundle.contents.count
    for i in range(count):
        entry = c_bundle.contents.entries[i]
        name = entry.name.decode("utf-8", errors="replace")
        mtype = ModelType(entry.type)
        m_inst = _wrap_c_model(entry.model, mtype)
        bundle.add_model(name, m_inst)

    # Free the C bundle container without freeing the model payloads (handed off to Python)
    lsci.DelModelBundle(ctypes.byref(c_bundle), 0)
    return bundle


def import_auto(filepath: str) -> Any | ModelBundle:
    """
    Automatically detects whether a file is a single model (.lsm) or a multi-model bundle (.lsb),
    and imports it appropriately.

    Returns:
        - For single model files: The model instance (e.g. PCA, PLS).
        - For multi-model bundle files: A ModelBundle instance containing all models.
    """
    info = inspect_model_file(filepath)
    if info["is_single"]:
        return import_model(filepath, expected_type=info["primary_type"])
    else:
        return import_model_bundle(filepath)
