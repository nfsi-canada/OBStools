project('rfpy', 'c',
  version : '0.2.0',
  license: 'MIT',
  meson_version: '>=0.64.0',
)

py_mod = import('python')
py = py_mod.find_installation()
py_dep = py.dependency()

py.install_sources(
  'obstools/__init__.py',
  'obstools/atacr/classes.py',
  'obstools/atacr/plotting.py',
  'obstools/atacr/utils.py',
  'obstools/atacr/__init__.py',
  'obstools/comply/classes.py',
  'obstools/comply/__init__.py',
  'obstools/scripts/atacr_download_data.py',
  'obstools/scripts/atacr_daily_spectra.py',
  'obstools/scripts/atacr_clean_spectra.py',
  'obstools/scripts/atacr_transfer_functions.py',
  'obstools/scripts/atacr_download_event.py',
  'obstools/scripts/atacr_correct_event.py',
  'obstools/scripts/comply_calculate.py',
  'obstools/scripts/__init__.py',
  pure: false,
  preserve_path: true,
  )
  