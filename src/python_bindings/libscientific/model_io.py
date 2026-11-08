# model_io.py libscientific python binding
#
# Copyright (C) 2026 designed, written and maintained by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
"""Saving and loading models in the libscientific binary formats.

Two file formats, both without external dependencies:

- ``.lsm``: a single model;
- ``.lsb``: a bundle of named models (e.g. a PCA and a PLS of the same
  project) in one file.

Most of the time the ``save``/``load`` methods of the models are enough
(``PCA.save("model.lsm")``). This module adds bundles and automatic
detection of the file content:

    >>> import os, tempfile
    >>> from libscientific.pca import PCA
    >>> from libscientific.pls import PLS
    >>> from libscientific.model_io import ModelBundle, import_auto
    >>> X = [[1.0, 2.0], [2.0, 1.0], [3.0, 4.0], [4.0, 3.0]]
    >>> pca = PCA(npc=2); pca.fit(X)
    >>> pls = PLS(nlv=2); pls.fit(X, [[1.0], [2.0], [3.0], [4.0]])
    >>> path = os.path.join(tempfile.mkdtemp(), "project.lsb")
    >>> ModelBundle({"pca": pca, "pls": pls}).save(path)
    >>> bundle = import_auto(path)
    >>> sorted(bundle.keys()), bundle["pca"].npc
    (['pca', 'pls'], 2)
"""

import ctypes
from enum import IntEnum
from typing import Any

from libscientific.loadlibrary import load_libscientific_library

lsci = load_libscientific_library()


class ModelType(IntEnum):
    """Type of a model stored in a model file.

    Members: ``UNKNOWN`` (0), ``PCA`` (1), ``PLS`` (2), ``CPCA`` (3), ``MLR`` (4),
    ``LDA`` (5), ``ICA`` (6), ``UPCA`` (7), ``UPLS`` (8), ``ANY`` (255, accept any
    type when importing).
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
    """Raise the Python exception matching a negative C status code."""
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
    """C structure filled by ``InspectModelFile`` (low-level API, see :func:`inspect_model_file`)."""

    _fields_ = [
        ("num_models", ctypes.c_uint32),
        ("is_single", ctypes.c_bool),
        ("primary_type", ctypes.c_uint32),
        ("primary_name", ctypes.c_char * 64),
        ("format_version", ctypes.c_uint32),
        ("total_size", ctypes.c_uint64),
    ]


class MODEL_ENTRY(ctypes.Structure):
    """C structure of one named model inside a bundle (low-level API)."""

    _fields_ = [
        ("name", ctypes.c_char * 64),
        ("type", ctypes.c_uint32),
        ("model", ctypes.c_void_p),
    ]


class MODEL_BUNDLE(ctypes.Structure):
    """C structure of a model bundle (low-level API, see :class:`ModelBundle`)."""

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


_CLASS_TO_TYPE = {
    "PCA": "PCA",
    "PLS": "PLS",
    "CPCA": "CPCA",
    "MLR": "MLR",
    "LDA": "LDA",
    "UPCA": "UPCA",
    "UPLS": "UPLS",
}


def _extract_c_model_and_type(
    obj: Any, explicit_type: ModelType | None = None
) -> tuple[ctypes.c_void_p, ModelType]:
    """
    Return the C pointer and :class:`ModelType` of a model.

    Parameters
    ----------
    obj : model object or ctypes pointer
        A high-level model (``PCA``, ``PLS``, ``CPCA``, ``MLR``, ``LDA``,
        ``UPCA``, ``UPLS``) or a pointer to the matching C structure
        (``PCAMODEL``, ...).
    explicit_type : ModelType, optional
        Type to use when it cannot be detected (e.g. a ``c_void_p``).

    Returns
    -------
    tuple[ctypes.c_void_p, ModelType]

    Raises
    ------
    ValueError
        If the type cannot be determined.
    """
    cls_name = type(obj).__name__
    if hasattr(obj, "model") and cls_name in _CLASS_TO_TYPE:
        return ctypes.cast(obj.model, ctypes.c_void_p), ModelType[_CLASS_TO_TYPE[cls_name]]
    if hasattr(obj, "_type_"):
        inner_name = obj._type_.__name__
        if inner_name.endswith("MODEL") and inner_name[:-5] in _CLASS_TO_TYPE:
            return ctypes.cast(obj, ctypes.c_void_p), ModelType[inner_name[:-5]]
    if explicit_type is not None and explicit_type != ModelType.UNKNOWN:
        return ctypes.cast(obj, ctypes.c_void_p), explicit_type
    raise ValueError(f"Unable to automatically determine model type for {obj} ({cls_name})")


def _wrap_c_model(ptr_val: Any, model_type: ModelType) -> Any:
    """
    Wrap a C model pointer, allocated by the library, into its Python class.

    The returned object owns the pointer and frees it when garbage collected.
    Preprocessing options that are not stored in the file (``scaling``,
    ``xscaling``, ``yscaling``) get their default values; they do not affect
    ``predict``, which uses the stored column averages and scalings.

    Parameters
    ----------
    ptr_val : int or ctypes pointer
        Address of the C model.
    model_type : ModelType
        Type of the model.

    Returns
    -------
    PCA, PLS, CPCA, MLR, LDA, UPCA, UPLS, or the raw pointer for other types
    """
    if not ptr_val:
        return None
    if model_type == ModelType.PCA:
        from libscientific.pca import PCA, PCAMODEL

        pca_ptr = ctypes.cast(ptr_val, ctypes.POINTER(PCAMODEL))
        pca_inst = PCA.__new__(PCA)
        pca_inst.model = pca_ptr
        pca_inst.scaling = 1
        pca_inst.npc = pca_ptr.contents.scores[0].col if bool(pca_ptr.contents.scores) else 0
        return pca_inst
    if model_type == ModelType.PLS:
        from libscientific.pls import PLS, PLSMODEL

        pls_ptr = ctypes.cast(ptr_val, ctypes.POINTER(PLSMODEL))
        pls_inst = PLS.__new__(PLS)
        pls_inst.model = pls_ptr
        pls_inst.xscaling = 1
        pls_inst.yscaling = 0
        pls_inst.nlv = pls_ptr.contents.xscores[0].col if bool(pls_ptr.contents.xscores) else 0
        return pls_inst
    if model_type == ModelType.CPCA:
        from libscientific.cpca import CPCA, CPCAMODEL

        cpca_ptr = ctypes.cast(ptr_val, ctypes.POINTER(CPCAMODEL))
        cpca_inst = CPCA.__new__(CPCA)
        cpca_inst.model = cpca_ptr
        cpca_inst.scaling = 1
        cpca_inst.npc = (
            cpca_ptr.contents.super_scores[0].col if bool(cpca_ptr.contents.super_scores) else 0
        )
        return cpca_inst
    if model_type == ModelType.MLR:
        from libscientific.mlr import MLR, MLRMODEL

        mlr_inst = MLR.__new__(MLR)
        mlr_inst.model = ctypes.cast(ptr_val, ctypes.POINTER(MLRMODEL))
        return mlr_inst
    if model_type == ModelType.LDA:
        from libscientific.lda import LDA, LDAMODEL

        lda_inst = LDA.__new__(LDA)
        lda_inst.model = ctypes.cast(ptr_val, ctypes.POINTER(LDAMODEL))
        return lda_inst
    if model_type == ModelType.UPCA:
        from libscientific.upca import UPCA, UPCAMODEL

        upca_ptr = ctypes.cast(ptr_val, ctypes.POINTER(UPCAMODEL))
        upca_inst = UPCA.__new__(UPCA)
        upca_inst.model = upca_ptr
        upca_inst.scaling = 1
        upca_inst.npc = upca_ptr.contents.scores[0].col if bool(upca_ptr.contents.scores) else 0
        return upca_inst
    if model_type == ModelType.UPLS:
        from libscientific.upls import UPLS, UPLSMODEL

        upls_ptr = ctypes.cast(ptr_val, ctypes.POINTER(UPLSMODEL))
        upls_inst = UPLS.__new__(UPLS)
        upls_inst.model = upls_ptr
        upls_inst.xscaling = 1
        upls_inst.yscaling = 0
        upls_inst.nlv = upls_ptr.contents.xscores[0].col if bool(upls_ptr.contents.xscores) else 0
        return upls_inst
    return ptr_val


class ModelBundle:
    """A dictionary-like collection of named models, saved to one ``.lsb`` file.

    Parameters
    ----------
    models : dict[str, model], list of models or None
        Initial models. A list gets the names ``model_0``, ``model_1``, ...

    Notes
    -----
    Supports ``bundle[name]``, ``bundle[name] = model``, ``del bundle[name]``,
    ``name in bundle``, ``len(bundle)``, iteration over the names, ``keys()``,
    ``values()`` and ``items()``.

    Examples
    --------
    >>> from libscientific.pca import PCA
    >>> pca = PCA(npc=1); pca.fit([[1.0, 2.0], [2.0, 1.0], [3.0, 4.0]])
    >>> bundle = ModelBundle({"pca": pca})
    >>> len(bundle), "pca" in bundle
    (1, True)
    """

    def __init__(self, models: dict[str, Any] | list[Any] | tuple[Any, ...] | None = None):
        """Create the bundle (see the class documentation)."""
        self._models: dict[str, Any] = {}
        if models:
            if isinstance(models, dict):
                for k, v in models.items():
                    self.add_model(k, v)
            elif isinstance(models, (list, tuple)):
                for idx, v in enumerate(models):
                    self.add_model(f"model_{idx}", v)

    def add_model(self, name: str, model: Any, model_type: ModelType | None = None):
        """Add (or replace) a named model.

        Parameters
        ----------
        name : str
            Name of the model inside the bundle.
        model : PCA, PLS, CPCA, MLR, LDA, UPCA, UPLS or C model pointer
            The model.
        model_type : ModelType, optional
            Needed only for raw C pointers whose type cannot be detected.
        """
        if not name or not isinstance(name, str):
            raise ValueError("Model name must be a non-empty string")
        self._models[name] = model

    def get_model(self, name: str) -> Any:
        """Return a model by name.

        Parameters
        ----------
        name : str
            Name of the model.

        Returns
        -------
        model

        Raises
        ------
        KeyError
            If there is no model with that name.
        """
        return self._models[name]

    def __getitem__(self, name: str) -> Any:
        """Return a model by name (``bundle[name]``)."""
        return self._models[name]

    def __setitem__(self, name: str, model: Any):
        """Add or replace a model (``bundle[name] = model``)."""
        self.add_model(name, model)

    def __delitem__(self, name: str):
        """Remove a model (``del bundle[name]``)."""
        del self._models[name]

    def __contains__(self, name: str) -> bool:
        """Return ``True`` if a model with that name exists."""
        return name in self._models

    def __len__(self) -> int:
        """Return the number of models."""
        return len(self._models)

    def __iter__(self):
        """Iterate over the model names."""
        return iter(self._models)

    def keys(self):
        """Return the model names."""
        return self._models.keys()

    def values(self):
        """Return the models."""
        return self._models.values()

    def items(self):
        """Return ``(name, model)`` pairs."""
        return self._models.items()

    def save(self, filepath: str):
        """Save all the models to a bundle file.

        Parameters
        ----------
        filepath : str
            Destination path (``.lsb`` recommended).

        Raises
        ------
        OSError, ValueError
            If a model cannot be stored or the file cannot be written.
        """
        save_model_bundle(filepath, self)

    def __repr__(self) -> str:
        names = ", ".join(repr(k) for k in self._models.keys())
        return f"<ModelBundle count={len(self._models)} models=[{names}]>"


def inspect_model_file(filepath: str) -> dict[str, Any]:
    """Read the header of a model file without loading the models.

    Parameters
    ----------
    filepath : str
        Path of a ``.lsm`` or ``.lsb`` file.

    Returns
    -------
    dict
        ``num_models`` (int), ``is_single`` (bool), ``primary_type``
        (:class:`ModelType` of the first model), ``primary_name`` (str),
        ``format_version`` (int) and ``total_size`` (bytes).

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    OSError
        If it is not a libscientific model file.
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
    """Return ``True`` if a file contains more than one model.

    Parameters
    ----------
    filepath : str
        Path of a ``.lsm`` or ``.lsb`` file.

    Returns
    -------
    bool
    """
    info = inspect_model_file(filepath)
    return not info["is_single"]


def save_model(
    filepath: str, model: Any, model_name: str = "default", model_type: ModelType | None = None
):
    """Save one model to a ``.lsm`` file.

    Equivalent to the ``save`` method of the models, with a custom name.

    Parameters
    ----------
    filepath : str
        Destination path.
    model : PCA, PLS, CPCA, MLR, LDA, UPCA, UPLS or C model pointer
        The model to save.
    model_name : str, default "default"
        Name stored in the file.
    model_type : ModelType, optional
        Needed only for raw C pointers whose type cannot be detected.

    Raises
    ------
    OSError, ValueError
        If the model cannot be saved.
    """
    c_ptr, mtype = _extract_c_model_and_type(model, model_type)
    fp = filepath.encode("utf-8")
    nm = model_name.encode("utf-8")
    st = lsci.SaveModel(fp, nm, int(mtype), c_ptr)
    _check_status(st, f"SaveModel('{filepath}')")


def import_model(
    filepath: str, expected_type: ModelType = ModelType.ANY, target_model: Any = None
) -> Any:
    """Load a single model from a file.

    Parameters
    ----------
    filepath : str
        Path of a ``.lsm`` file (or a ``.lsb`` file with one model).
    expected_type : ModelType, default ModelType.ANY
        Raise ``ValueError`` if the file contains another type of model.
    target_model : model, optional
        An existing model object to fill; by default a new object of the right
        class is created.

    Returns
    -------
    PCA, PLS, CPCA, MLR, LDA, UPCA or UPLS
        The loaded model.

    Raises
    ------
    ValueError
        If the type does not match ``expected_type``, or the file contains
        several models (use :func:`import_model_bundle` or :func:`import_auto`).
    FileNotFoundError, OSError
        If the file cannot be read.
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
            if cpca_inst.model and cpca_inst.model.contents.super_scores:
                cpca_inst.npc = cpca_inst.model.contents.super_scores[0].col
            return cpca_inst
        else:
            out_model = ctypes.c_void_p(None)
            out_type = ctypes.c_uint32(0)
            out_bundle = ctypes.POINTER(MODEL_BUNDLE)()
            st = lsci.ImportModelAuto(
                fp, ctypes.byref(out_model), ctypes.byref(out_type), ctypes.byref(out_bundle)
            )
            # ImportModelAuto returns the number of models read (>= 1) or a
            # negative error code.
            if st < 0:
                _check_status(st, f"ImportModelAuto('{filepath}')")
            if out_model.value:
                return _wrap_c_model(out_model.value, ModelType(out_type.value))
            if out_bundle:
                lsci.DelModelBundle(ctypes.byref(out_bundle), 1)
            raise ValueError(
                f"'{filepath}' contains {st} models: use import_model_bundle() or import_auto()"
            )


def save_model_bundle(
    filepath: str, bundle: ModelBundle | dict[str, Any] | list[Any] | tuple[Any, ...]
):
    """Save several models to one ``.lsb`` file.

    Parameters
    ----------
    filepath : str
        Destination path.
    bundle : ModelBundle, dict[str, model] or list of models
        The models; a list gets the names ``model_0``, ``model_1``, ...

    Raises
    ------
    TypeError
        If ``bundle`` has an unsupported type.
    OSError, ValueError
        If a model cannot be stored or the file cannot be written.
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
    """Load all the models of a bundle file.

    Parameters
    ----------
    filepath : str
        Path of a ``.lsb`` file.

    Returns
    -------
    ModelBundle
        The models, by name.

    Raises
    ------
    OSError
        If the file cannot be read.
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
    """Load a model file, whatever it contains.

    Parameters
    ----------
    filepath : str
        Path of a ``.lsm`` or ``.lsb`` file.

    Returns
    -------
    model or ModelBundle
        The model for single-model files, a :class:`ModelBundle` otherwise.
    """
    info = inspect_model_file(filepath)
    if info["is_single"]:
        return import_model(filepath, expected_type=info["primary_type"])
    else:
        return import_model_bundle(filepath)
