.. _contributing:

*********************
Contributing to pyque
*********************

.. contents:: Table of contents:
   :local:

The ``meta`` subpackage is very important and is the basic of our
model.

-  The ``param`` module defines a ``Parameter`` class that is the
   parameter of the namelist, e.g., for ``pw.x``, ``calculation`` is a
   parameter, and it can takes 7 different values (``'scf'``,
   ``'nscf'``, etc.). When we read a PWscf input file, a string is throw
   to us—any one of the seven. Well, becuase this parameter has a type
   ``'str'``, the string, supposed to be ``'scf'``, is converted to be a
   string, nothing change, of course. But if we read parameter
   ``ibrav``, then the type will be an ``'int'``, so the input, like
   ``'0'``, will be converted to ``0``.
-  The ``tree`` module defines a data structure that will be used as a
   superclass for our ``input`` object.
-  In ``pwscf`` I defined an object ``SCFStandardInput``, which is a
   standard file that can be read directly by Quantum ESPRESSO. It needs
   to be built by ``SCFInputReader,build_input_tree`` method. It also
   has a ``to_text_file`` method, because I want it can be rewrite back
   to a file once we have modified a parameter in it.

Please read ``PhononInputReader`` in ``readers.phonon`` and
``basics.phonon_params`` and ``NamelistReader`` in ``readers.simple`` to
get understood what I am talking about. ‘namelist’ is the same as ‘card’
and used interchangablly, but I will correct all of them to be
‘namelist’ later, since this might be defined by Quantum ESPRESSO.

The ``miscellaneous`` subpackage should store something very useful but
not crucial to our tasks.

The ``shell`` subpackage store pure shell scripts and should be
rewritten to Python, then discarded.

The ``plotters`` and ``readers`` subpackages defined a lot to plot and
read files.

The ``calculator`` subpackage is used for QE-related calculation.

Every subpackage (what functions) should have a ``__init__.py`` file,
but not the unit test folders.
