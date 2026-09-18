##############################################################################
## Copyright (c) 2017 Jeet Sukumaran.
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.
##     * Redistributions in binary form must reproduce the above copyright
##       notice, this list of conditions and the following disclaimer in the
##       documentation and/or other materials provided with the distribution.
##     * The names of its contributors may not be used to endorse or promote
##       products derived from this software without specific prior written
##       permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
## IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
## THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
## PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL JEET SUKUMARAN BE LIABLE FOR ANY
## DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
## (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
## LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
## AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
## SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
##############################################################################

from __future__ import annotations

import codecs
import collections
import collections.abc
import csv
import logging
import os
import re
import shutil
import sys
import tempfile
from io import StringIO  # noqa: F401  (оставлен для совместимости импортов)

__all__ = [
    "open_output_file_for_csv_writer",
    "get_csv_writer",
    "write_dict_csv",
    "parse_legacy_configuration",
    "bytes_to_text",
    "communicate_process",
    "parse_fieldname_and_value",
    "filter_columns_using_master_template_file",
    "extract_fieldnames_from_file",
    "filter_columns_from_file",
    "RunLogger",
    "CaseInsensitiveDict",
]


##############################################################################
## CSV File Handling

def open_output_file_for_csv_writer(filepath, is_append: bool = False):

    if filepath is None or filepath == "-":
        return sys.stdout, False
    mode = "a" if is_append else "w"
    stream = open(filepath, mode, newline="", encoding="utf-8")
    return stream, True


# Старое имя оставлено как алиас — на случай, если кто-то ещё его вызывает.
open_destput_file_for_csv_writer = open_output_file_for_csv_writer


def get_csv_writer(
        dest,
        fieldnames=None,
        delimiter: str = "\t",
        lineterminator: str | None = None,
        ):
   
    if lineterminator is None:
        lineterminator = os.linesep
    if isinstance(dest, (str, bytes, os.PathLike)):
        dest, _ = open_output_file_for_csv_writer(filepath=dest)
    return csv.DictWriter(
            dest,
            fieldnames=fieldnames,
            restval="NA",
            delimiter=delimiter,
            lineterminator=lineterminator,
            )


def write_dict_csv(
        list_of_dicts,
        filepath,
        fieldnames=None,
        is_no_header_row: bool = False,
        is_append: bool = False,
        ):
 
    if not list_of_dicts:
        return
    if fieldnames is None:
        fieldnames = list(list_of_dicts[0].keys())

    dest, should_close = open_output_file_for_csv_writer(
            filepath=filepath,
            is_append=is_append,
            )
    try:
        writer = get_csv_writer(
                dest=dest,
                fieldnames=fieldnames,
                delimiter=",",
                )
        if not is_no_header_row:
            writer.writeheader()
        writer.writerows(list_of_dicts)
    finally:
        if should_close:
            dest.close()


##############################################################################
## Configuration File Handling

def parse_legacy_configuration(filepath, config_d=None):
    recognized_preamble_keys = collections.OrderedDict([
            ("concentrationshape", float),
            ("concentrationscale", float),
            ("thetashape", float),
            ("thetascale", float),
            ("ancestralthetashape", float),
            ("ancestralthetascale", float),
            ("thetaparameters", str),
            ("taushape", float),
            ("tauscale", float),
            ("timeinsubspersite", float),
            ("bottleproportionshapea", float),
            ("bottleproportionshapeb", float),
            ("bottleproportionshared", float),
            ("migrationshape", float),
            ("migrationscale", float),
            ("numtauclasses", float),
            ])
    sample_table_keys = [
            ("taxon_label", str),
            ("locus_label", str),
            ("ploidy_factor", float),
            ("mutation_rate_factor", float),
            ("num_genes_deme0", int),
            ("num_genes_deme1", int),
            ("ti_tv_rate_ratio", float),
            ("num_sites", int),
            ("freq_a", float),
            ("freq_c", float),
            ("freq_g", float),
            ("alignment_filepath", str),
            ]
    sample_table_begin_pattern = re.compile(r"^\s*BEGIN\s+SAMPLE_TBL\s*$", re.I)
    sample_table_end_pattern = re.compile(r"^\s*END\s+SAMPLE_TBL\s*$", re.I)
    sample_table_splitter = re.compile(r"\s+", re.I)

    if config_d is None:
        config_d = CaseInsensitiveDict()
    config_d["params"] = {}
    config_d["locus_info"] = []
    section = "preamble"

    with open(filepath, encoding="utf-8") as src:
        for row_idx, row in enumerate(src):
            row = row.strip()
            if not row:
                continue
            comment_start_idx = row.find("#")
            if section == "sample-table" and comment_start_idx > 0:
                raise ValueError(
                    "Configuration file '{}', row {}: sample table section "
                    "cannot contain mid-row comment".format(filepath, row_idx + 1)
                )
            if comment_start_idx > -1:
                row = row[0:comment_start_idx]
            row = row.strip()
            if not row:
                continue
            if section == "preamble":
                if sample_table_begin_pattern.match(row):
                    section = "sample-table"
                    continue
                row_parts = row.split("=")
                if len(row_parts) != 2:
                    raise ValueError(
                        "Configuration file '{}', row {}: expected 'key = value'".format(
                            filepath, row_idx + 1
                        )
                    )
                key = row_parts[0].strip()
                case_normalized_key = key.lower()
                if case_normalized_key not in recognized_preamble_keys:
                    raise ValueError(
                        "Configuration file '{}', row {}: unrecognized preamble key '{}'".format(
                            filepath, row_idx + 1, key
                        )
                    )
                config_d["params"][key] = recognized_preamble_keys[case_normalized_key](
                        row_parts[1].strip()
                )
            else:
                if sample_table_end_pattern.match(row):
                    continue
                cols = sample_table_splitter.split(row)
                if len(cols) != len(sample_table_keys):
                    raise ValueError(
                        "Configuration file '{}', row {}: expecting {} columns "
                        "but only found {}".format(
                            filepath, row_idx + 1, len(sample_table_keys), len(cols)
                        )
                    )
                locus_info = {}
                for (key, val_type), val in zip(sample_table_keys, cols):
                    locus_info[key] = val_type(val)
                config_d["locus_info"].append(locus_info)
    return config_d


##############################################################################
## Process Control/Handling

try:
    ENCODING = sys.getdefaultencoding() or "utf-8"
except Exception:  # pragma: no cover
    ENCODING = "utf-8"

if not ENCODING:
    ENCODING = "utf-8"


def bytes_to_text(s):
    if s is None:
        return None
    if isinstance(s, str):
        return s
    return codecs.decode(s, ENCODING)


def communicate_process(p, commands=None, timeout=None):
    if isinstance(commands, (list, tuple)):
        commands = "\n".join(str(c) for c in commands)
    if commands is not None:
        commands = commands.encode(ENCODING)
    if timeout is None:
        stdout, stderr = p.communicate(commands)
    else:
        stdout, stderr = p.communicate(commands, timeout=timeout)
    return bytes_to_text(stdout), bytes_to_text(stderr)


##############################################################################
## Command line processing

def parse_fieldname_and_value(labels):
    if not labels:
        return collections.OrderedDict()
    fieldname_value_map = collections.OrderedDict()
    for label in labels:
        match = re.match(r"\s*(.*?)\s*:\s*(.*)\s*", label)
        if not match:
            raise ValueError(
                "Cannot parse fieldname and label "
                "(format required: fieldname:value): {}".format(label)
            )
        fieldname, value = match.groups(0)
        fieldname_value_map[fieldname] = value
    return fieldname_value_map


##############################################################################
## Post-processing

def extract_fieldnames_from_file(src, field_delimiter: str = "\t"):
    """Возвращает список имён колонок из файла или открытого потока."""
    if isinstance(src, (str, bytes, os.PathLike)):
        with open(src, encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f, delimiter=field_delimiter, quoting=csv.QUOTE_NONE)
            return list(reader.fieldnames or [])
    reader = csv.DictReader(src, delimiter=field_delimiter, quoting=csv.QUOTE_NONE)
    return list(reader.fieldnames or [])


def filter_columns_from_file(
        src,
        dest,
        columns_to_retain,
        field_delimiter: str = "\t",
        ):
    src_is_path = isinstance(src, (str, bytes, os.PathLike))
    dest_is_path = isinstance(dest, (str, bytes, os.PathLike))

    src_stream = open(src, encoding="utf-8", newline="") if src_is_path else src
    dest_stream = open(dest, "w", encoding="utf-8", newline="") if dest_is_path else dest

    try:
        source_reader = csv.DictReader(
                src_stream,
                delimiter=field_delimiter,
                quoting=csv.QUOTE_NONE,
                )
        target_writer = csv.DictWriter(
                dest_stream,
                delimiter=field_delimiter,
                quoting=csv.QUOTE_NONE,
                fieldnames=[],
                restval="NA",
                lineterminator=os.linesep,
                )
        to_delete = None
        to_keep = []
        for row in source_reader:
            if to_delete is None:
                source_fields = source_reader.fieldnames or []
                to_delete = {k for k in source_fields if k not in columns_to_retain}
                to_keep = [k for k in columns_to_retain if k in source_fields]
                target_writer.fieldnames = to_keep
                target_writer.writeheader()
            for key in to_delete:
                row.pop(key, None)
            target_writer.writerow(row)
    finally:
        if src_is_path:
            src_stream.close()
        if dest_is_path:
            dest_stream.close()


def filter_columns_using_master_template_file(
        dest,
        master_file,
        source_file,
        field_delimiter: str = "\t",
        ):
    master_field_names = extract_fieldnames_from_file(
            src=master_file, field_delimiter=field_delimiter
            )
    filter_columns_from_file(
            dest=dest,
            src=source_file,
            columns_to_retain=master_field_names,
            field_delimiter=field_delimiter,
            )


##############################################################################
## Logging

_LOGGING_LEVEL_ENVAR = "GERENUK_LOGGING_LEVEL"
_LOGGING_FORMAT_ENVAR = "GERENUK_LOGGING_FORMAT"


class RunLogger(object):

    NOTSET_MESSAGING_LEVEL = logging.NOTSET
    DEBUG_MESSAGING_LEVEL = logging.DEBUG
    INFO_MESSAGING_LEVEL = logging.INFO
    WARNING_MESSAGING_LEVEL = logging.WARNING
    ERROR_MESSAGING_LEVEL = logging.ERROR
    CRITICAL_MESSAGING_LEVEL = logging.CRITICAL

    def __init__(self, **kwargs):
        self.name = kwargs.get("name", "RunLog")
        self._log = logging.getLogger(self.name)
        self._log.setLevel(RunLogger.DEBUG_MESSAGING_LEVEL)
        self.handlers = []

        if kwargs.get("log_to_stderr", True):
            handler1 = logging.StreamHandler()
            stderr_logging_level = self.get_logging_level(
                    kwargs.get("stderr_logging_level", RunLogger.INFO_MESSAGING_LEVEL)
            )
            handler1.setLevel(stderr_logging_level)
            handler1.setFormatter(self.get_default_formatter())
            self._log.addHandler(handler1)
            self.handlers.append(handler1)

        if kwargs.get("log_to_file", True):
            if "log_stream" in kwargs:
                log_stream = kwargs.get("log_stream")
            else:
                log_stream = open(
                        kwargs.get("log_path", self.name + ".log"),
                        "w",
                        encoding="utf-8",
                        )
            handler2 = logging.StreamHandler(log_stream)
            file_logging_level = self.get_logging_level(
                    kwargs.get("file_logging_level", RunLogger.DEBUG_MESSAGING_LEVEL)
            )
            handler2.setLevel(file_logging_level)
            handler2.setFormatter(self.get_default_formatter())
            self._log.addHandler(handler2)
            self.handlers.append(handler2)

        self._system = None

    # --- system property ---------------------------------------------------

    def _get_system(self):
        return self._system

    def _set_system(self, system):
        self._system = system
        if self._system is None:
            formatter = self.get_default_formatter()
        else:
            formatter = self.get_simulation_generation_formatter()
        for handler in self.handlers:
            handler.setFormatter(formatter)

    system = property(_get_system, _set_system)

    # --- levels ------------------------------------------------------------

    def get_logging_level(self, level=None):
        levels = {
            "NOTSET": RunLogger.NOTSET_MESSAGING_LEVEL,
            "DEBUG": RunLogger.DEBUG_MESSAGING_LEVEL,
            "INFO": RunLogger.INFO_MESSAGING_LEVEL,
            "WARNING": RunLogger.WARNING_MESSAGING_LEVEL,
            "ERROR": RunLogger.ERROR_MESSAGING_LEVEL,
            "CRITICAL": RunLogger.CRITICAL_MESSAGING_LEVEL,
        }
        if isinstance(level, int) and level in levels.values():
            return level
        if level is not None:
            level_name = str(level).upper()
        elif _LOGGING_LEVEL_ENVAR in os.environ:
            level_name = os.environ[_LOGGING_LEVEL_ENVAR].upper()
        else:
            level_name = "NOTSET"
        return levels.get(level_name, RunLogger.NOTSET_MESSAGING_LEVEL)

    # --- formatters --------------------------------------------------------

    def get_default_formatter(self):
        f = logging.Formatter("[%(asctime)s] %(message)s")
        f.datefmt = "%Y-%m-%d %H:%M:%S"
        return f

    def get_simulation_generation_formatter(self):
        f = logging.Formatter("[%(asctime)s] %(message)s")
        f.datefmt = "%Y-%m-%d %H:%M:%S"
        return f

    def get_rich_formatter(self):
        f = logging.Formatter("[%(asctime)s] [%(levelname)s] %(message)s")
        f.datefmt = "%Y-%m-%d %H:%M:%S"
        return f

    def get_simple_formatter(self):
        f = logging.Formatter("%(message)s")
        return f

    def get_raw_formatter(self):
        return logging.Formatter("%(message)s")

    def get_logging_formatter(self, fmt=None):
        if fmt is not None:
            fmt = fmt.upper()
        elif _LOGGING_FORMAT_ENVAR in os.environ:
            fmt = os.environ[_LOGGING_FORMAT_ENVAR].upper()
        else:
            fmt = "DEFAULT"

        if fmt == "RICH":
            formatter = self.get_rich_formatter()
        elif fmt == "SIMPLE":
            formatter = self.get_simple_formatter()
        elif fmt == "NONE":
            formatter = self.get_raw_formatter()
        else:
            formatter = self.get_default_formatter()

        formatter.datefmt = "%H:%M:%S"
        return formatter

    # --- logging helpers ---------------------------------------------------

    def debug(self, msg):
        self._log.debug("[DEBUG] %s", msg)

    def info(self, msg):
        self._log.info(msg)

    def warning(self, msg):
        self._log.warning(msg)

    def error(self, msg):
        self._log.error(msg)

    def critical(self, msg):
        self._log.critical(msg)

    def log(self, msg, level):
        if level == RunLogger.DEBUG_MESSAGING_LEVEL:
            self.debug(msg)
        elif level == RunLogger.INFO_MESSAGING_LEVEL:
            self.info(msg)
        elif level == RunLogger.WARNING_MESSAGING_LEVEL:
            self.warning(msg)
        elif level == RunLogger.ERROR_MESSAGING_LEVEL:
            self.error(msg)
        elif level == RunLogger.CRITICAL_MESSAGING_LEVEL:
            self.critical(msg)
        else:
            raise ValueError("Unrecognized messaging level: {}".format(level))

    def flush(self):
        for handler in self.handlers:
            handler.flush()


###############################################################################
# CaseInsensitiveDict
#
# From: https://github.com/kennethreitz/requests
# Copyright 2014 Kenneth Reitz, Apache License 2.0
###############################################################################

class CaseInsensitiveDict(collections.abc.MutableMapping):

    def __init__(self, data=None, **kwargs):
        self._store = {}
        if data is None:
            data = {}
        self.update(data, **kwargs)

    def __setitem__(self, key, value):
        self._store[key.lower()] = (key, value)

    def __getitem__(self, key):
        return self._store[key.lower()][1]

    def __delitem__(self, key):
        del self._store[key.lower()]

    def __iter__(self):
        return (casedkey for casedkey, _ in self._store.values())

    def __len__(self):
        return len(self._store)

    def lower_items(self):
        return ((lowerkey, keyval[1]) for lowerkey, keyval in self._store.items())

    def __eq__(self, other):
        if isinstance(other, collections.abc.Mapping):
            other = CaseInsensitiveDict(other)
        else:
            return NotImplemented
        return dict(self.lower_items()) == dict(other.lower_items())

    def copy(self):
        return CaseInsensitiveDict(self._store.values())

    def __repr__(self):
        return "%s(%r)" % (self.__class__.__name__, dict(self.items()))
