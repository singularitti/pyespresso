# Type annotation

[TOC]

Python 3 does allow you to do type annotations. Simply like this:

```python
def func(arg1: str, arg2: float) -> float
	# do something
```

where the first argmument is `str` type, and second is `float` type, the type of returned value is `float`. These annotations does not have syntaxical meaning, they are just "annotations" for reading the code.

## List of lists

I adopt something like

```python
[[1, 2], [2, 3]]
[[float], [float]]
```

Note that here `[float]` does not mean there is only one `float` inside `list`, i.e., this does not show its "shape" but only type, for simplicity. Similarly, `[[float]]` does not mean there is only one `[float]` inside the upper level `list`.