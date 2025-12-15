from .equations import *
from .equation_container import EquationSet

EQUATION_REGISTRY = {
    "1": SimplestEquation,
    "2": HomogeneousEquation,
    "3": SumToProductEquation,
    "4": GroupingEquation,
    "5": PowerReductionEquation,
    "6": QuadraticTrigEquation,
    "7": DoubleAngleToQuadraticEquation,
    "8": LinearCombinationEquation,
    "9": ReducibleToHomogeneousEquation,
    "10": SymmetricEquation,
    "11": TanSubstitutionEquation,
    "12": SumTanCotanEquation,
    "13": BoundedSumEquation,
    "14": InverseTrigEquation,
}

__all__ = [
    'EquationSet',
    'EQUATION_REGISTRY',
]