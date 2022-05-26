Inheritance in Python
---------------------

``SeisFlows`` is built upon the object-oriented programming concept of
**inheritance**. This documentation page is a simple introduction to
this concept to help new users and developers understand how SeisFlows
is built and expected to operate.

**Inheritance** is the ability of one class to derive attributes from
another class, improving code *re-usability*. Some terminology used in
to talk about this inheritance is defined here:

-  **Base** (Baseclass): The foundational **Base** class which defines
   standard behavior. The Baseclass does not inherit any of its
   attributes or behavior.
-  **Parent** (Superclass): A class which is being inherited from. A
   Parent can be a Baseclass, but inheritance can also be daisy-chained.
-  **Child** (Subclass): A class that inherits some or all of its
   attributes from a parent.

Consider the following toy example where we define a **Base** class
which has some internal attributes and functions.

.. code:: ipython3

    class Base:
        """
        A Baseclass example. All SeisFlows modules contain a Base class which 
        defines the foundational structure that all inherited classes will adopt
        """
        def __init__(self, example_integer=5, example_float=1.2):
            """
            The init function defines instance-variables and their
            default values
            
            :type 
            """
            self.example_integer = example_integer
            self.example_float = example_float
            
        def check(self):
            """
            Check functions ensure that parameters are set correctly.
            This toy problem check function simply asserts types are
            set correctly.
            """
            assert(self.example_integer < 10), \
                "The example integer must be < 10"
    
            assert(self.example_float > 1.), \
                "The example float must be > 1."
            
        def manipulate(self):
            """
            Manipulate internal attributes
            
            Each module provides functions which serve a purpose in 
            the larger workflow. 
            
            :rtype: float
            :return: example integer added to example float
            """
            return self.example_integer + self.example_float

We can quickly look at the behavior of this class by creating an
instance of it.

.. code:: ipython3

    module = Base(example_integer=11, example_float=3.)
    module.check()


::


    ---------------------------------------------------------------------------

    AssertionError                            Traceback (most recent call last)

    <ipython-input-48-0d3a45ed832e> in <module>
          1 module = Base(example_integer=11, example_float=3.)
    ----> 2 module.check()
    

    <ipython-input-47-19eac124f82c> in check(self)
         21         """
         22         assert(self.example_integer < 10), \
    ---> 23             "The example integer must be < 10"
         24 
         25         assert(self.example_float > 1.), \


    AssertionError: The example integer must be < 10


.. code:: ipython3

    module = Base(example_integer=6, example_float=3.)
    module.check()
    print(module.manipulate())


.. parsed-literal::

    9.0


What problem does inheritance solve?
------------------------------------

Say we want to use the structure of the ``Baseclass``, but we need the
``manipulate`` function to return the subtraction of ``example_integer``
and ``example_float``, instead of their addition. There are a few ways
to approach this problem.

Plausible but cumbersome approaches
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  **Copy-paste**: One method of doing this would be to copy-paste the
   *entire* ``Baseclass`` (e.g., as ``BaseCopy``) and re-define the
   ``manipulate`` function. This would isntantly double our code length,
   with a lot of the new code being completly redundant. Additionally,
   if we made **any** changes to the ``Baseclass``, we would need to
   also make those changes to ``BaseCopy`` to keep functionality
   consistent.

-  **Create a new function**: Another method would be to define a
   completly new function, e.g., ``manipulate2``. This is more
   acceptable, BUT if some other script, function or module calls
   ``Base.manipulate()``, we will now need to make them call
   ``Base.manipulate2()`` instead. This involves a signficant amount of
   work. Similarly, consider the case where we want to go back to the
   original ``manipulate`` function.

The inheritance approach
~~~~~~~~~~~~~~~~~~~~~~~~

**Inheritance** solves this problem but allowing us to overwrite the
manipulate function by creating a ``Child`` class, which inherits the
properties of its ``Parent``. This results in the least amount of code
writing, keeps behavior consistent, and allows flexibility in editing
established code (e.g., the ``Baseclass``). Let’s see how this is done:

.. code:: ipython3

    class Super(Base):
        """
        This Superclass will now inherit all of the attributes of the Baseclass.
        It does nothing new.
        """
        pass

.. code:: ipython3

    module = Super(example_integer=6, example_float=3.)
    module.check()
    print(module.manipulate())


.. parsed-literal::

    9.0


Overwriting functions
~~~~~~~~~~~~~~~~~~~~~

To solve the problem stated above, we can totally overwrite the
manipulate function to provide different behavior

.. code:: ipython3

    class Super(Base):
        """
        This Superclass overwrites the manipulate function
        """
        def manipulate(self):
            """
            Manipulate internal attributes 
            
            :rtype: float
            :return: example integer subtracted from example float
            """ 
            return self.example_integer - self.example_float

.. code:: ipython3

    module = Super(example_integer=6, example_float=3.)
    module.check()
    print(module.manipulate())


.. parsed-literal::

    3.0


super() functions
~~~~~~~~~~~~~~~~~

The `super() <https://docs.python.org/3/library/functions.html#super>`__
function “returns a proxy object that delegates method calls to a parent
or sibling class.” In other words, super() calls the Parent class.

We can use the Python super() function to directly incorporate functions
from the parent class, allowing us to build upon previously written
code. This is useful if you don’t want to completely overwrite a
previously-defined function.

.. code:: ipython3

    class Super(Base):
        """
        This Superclass overwrites the manipulate function
        """
        def manipulate(self):
            """
            Manipulate internal attributes 
            
            :rtype: float
            :return: example integer subtracted from example float
            """ 
            added_values = super().manipulate()  # This calls Base.manipulate()
            print(f"added_values={added_values}")
            return added_values ** 2

.. code:: ipython3

    module = Super(example_integer=6, example_float=3.)
    module.check()
    print(module.manipulate())


.. parsed-literal::

    added_values=9.0
    81.0


Multiple inheritance
--------------------

Inheritance can be chained, meaning former ``Children`` can become
``Parents``! Although chaining inheritance can quickly become messy and
confusing, it is useful for extending existing capabilities without
having to make direct edits to the ``Parent`` classes.

Let’s say you want to inherit all of the capabilities of the Super
class, but you want to extend it further for your own specific workflow.
Here we define a ``Superer`` class, which inherits and extends the
``Super`` class (which itself inherits from the ``Base`` class).

.. code:: ipython3

    class Superer(Super):
        """
        This Superclass inherits from the Super class, which itself inherits from the Base class
        """
        def __init__(self, new_value=8, **kwargs):
            """
            We can extend the internal attributes in our Superclass. 
            The **kwargs allow us to be lazy and assume that the User understands class values must be
            passed all the way to the Baseclass
            """
            super().__init__(**kwargs)
            self.new_value = new_value
            
        def check(self):
            """
            We would like to extend the check function to address our new value,
            while still checking the original values
            """
            super().check()
            assert(self.new_value != 0), "New value must be > 0"
            
        def manipulate(self):
            """
            We can further manipulate this function, which itself has been changed in
            the Superclass.
            
            :rtype: float
            :return: example integer subtracted from example float
            """ 
            squared_values = super().manipulate()
            print(f"squared_values={squared_values}")
            return squared_values / 2
        
        def manipulate_more(self):
            """
            We can also define completely new functions which are not present in any of the Parent classes.
            This is useful when your Superclass needs to fully extend the functionalities of its Parents.
            """
            manipulated_value = self.manipulate()
            return self.new_value + manipulated_value

.. code:: ipython3

    module = Superer(example_integer=6, example_float=3., new_value=2)
    module.check()
    print(f"manipulate: {module.manipulate()}")
    print(f"manipulate_more: {module.manipulate_more()}")


.. parsed-literal::

    added_values=9.0
    squared_values=81.0
    manipulate: 40.5
    added_values=9.0
    squared_values=81.0
    manipulate_more: 42.5

