from inspect import isclass

from chemsmart.utils.mixins import RegistryMixin


class TestMixins:

    def test_get_subclasses(self):
        class Animal(RegistryMixin):
            pass

        class Mammal(Animal):
            pass

        class Reptile(Animal):
            pass

        class Dog(Mammal):
            pass

        subclasses = set(Animal.subclasses())
        # Expected subclasses
        expected = {Mammal, Reptile, Dog}

        # Test: Check if all expected subclasses are registered
        assert (
            subclasses == expected
        ), f"Expected {expected}, but got {subclasses}"

        # Test: Ensure all items in `subclasses()` are classes
        for cls in subclasses:
            assert isclass(cls), f"{cls} is not a class"

        # Test: Ensure all items are subclasses of Animal
        for cls in subclasses:
            assert issubclass(
                cls, Animal
            ), f"{cls} is not a subclass of Animal"
