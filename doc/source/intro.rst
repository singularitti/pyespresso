Intro to Data Structures
========================

In order to provide this package with greater extensibility, we build it on
a very clear hierachy. So users will not change much their code if we change
the implementations in the future.

``*Input`` objects
------------------

Quantum ESPRESSO input objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A simplest example is the Quantum ESPRESSO (QE) standard input files, here
we will parse them into ``*Input`` objects, so-called QE input objects.

A typical QE input object consists of two major parts:

1. namelists
2. cards

A typical namelist looks like this::

  &control
    calculation='scf'
    restart_mode='from_scratch',
    prefix='silicon',
    lelfield=.true.,
    nberrycyc=1
    pseudo_dir='$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
  /

Here I call each key-value pairs like ``calculation='scf'`` as a *parameter*,
this concept will be used later. The keys, like ``calculation`` are called
*name*s, and the values, like ``'scf'`` are called *value*s.
In short, namelists are built on top of parameters.

So, in the code, we define a very basic data structure, i.e.,
``NamelistVariable``. It has several attributes, including
``name``, like ``'calculation'``, which **must** be a string;
and ``default_type``, also a string, a valid value can be
``'str'``, ``'int'``, ``'float'``, ``'bool'``;
and ``value``, can be any of ``str``, ``int``, ``float``, ``bool``;
and ``default_value``, which means the default value given by Quantum ESPRESSO
documentation;
and ``in_namelist``, which shows what namelist this parameter should be in,
this is also a ``str`` type.

It also includes some methods, like
``to_qe_str``, which will convert its ``value`` attribute to a legal
Quantum ESPRESSO string for writing to a standard input file which can be
read by Quantum ESPRESSO. For example, a Python string will
be returned without change, a Python integer ``n`` will be convert to
``str(n)``, a Python float will be convert to a double precision Fortran
float, and a ``True`` or a ``False`` will be converted to ``'.true.'``
and ``'.false.'`` respectively.
And an instance of ``NamelistVariable`` can be also converted to a Python
``dict`` object by calling its ``to_dict()`` method, this is useful when
we want to save our data in new file format like JSON.

In a namelist of a QE input object, any parameter could be
directly accessed by dot access. For example,
we have an PWscf input object ``pw``, there is a
``'calculation'`` name in its ``control_namelist`` attribute,
the we can get it by

.. ipython::

  @verbatim
  In [1]: pw.control_namelist.calculation
  Out[1]: 'scf'

Simple and easy, huh? We can also set this by

.. ipython::

  @verbatim
  In [1]: pw.control_namelist.calculation = 'si'
  @verbatim
  

But remember, what you see may be just a value like ``'scf'``,
what actually you get is a ``CONTROLNamelistVariable``! We
shadow its implementation, but for power users, you should
know this point, and you can get any of the attributes we
introduced above for a ``NamelistVariable``!

For more information, refer the source code in ``pyque.meta.namelist``.

``BatchInput`` objects
^^^^^^^^^^^^^^^^^^^^^^

A ``BatchInput`` is mainly made of 5 parts:

1. shebang
2. directives
3. modules
4. commands
5. comments
