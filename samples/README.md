# Samples

Scripts with sample uses of the lib:

- **aligned.py**: aligned Mexico turbine, steady solution for radial forces distribution.
- **iea15mw.py**: similar to `aligned.py` but for IEA15MW reference turbine.
- **yaw.py**: Mexico turbine under yaw, three sections for multiple revolutions.
- **pitch_maneuver.py**: pitch maneauver, time evolution of forces for single
  section at fixed azimuth.

The results subfolder is ignored by git.

## How to

You need to include the python folder to use the sample scripts. This can be
done before calling the scripts or inside it:

```bash
export PYTHONPATH=/path/to/repo:$PYTHONPATH
```

OR

```python
import sys
sys.path.append('/path/to/repo')
```

If the Qt error is encoutered (`Qt: Session management error: None of the
authentication protocols specified are supported`), you can change the
matplotlib backend. For example:

```bash
MPLBACKEND='TkAgg' python pitch_maneuver.py
```