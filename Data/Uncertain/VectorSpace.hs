{-# LANGUAGE TypeFamilies                  #-}
{-# LANGUAGE ConstraintKinds               #-}


module Data.Uncertain.VectorSpace where


import Data.Basis
import Data.VectorSpace

import Data.Packed.Matrix
import Numeric.LinearAlgebra



class (HasBasis v) => FinDimVecSpace v where
  dimension :: v -> Int

instance FinDimVecSpace Double where
  dimension = const 1
instance FinDimVecSpace Float where
  dimension = const 1
instance ( FinDimVecSpace a, FinDimVecSpace b
         , s ~ Scalar a, s ~ Scalar b         )
         => FinDimVecSpace (a,b) where
  dimension (a,b) = dimension a + dimension b
instance ( FinDimVecSpace a, FinDimVecSpace b, FinDimVecSpace c
         , s ~ Scalar a, s ~ Scalar b, s ~ Scalar c             )
         => FinDimVecSpace (a,b,c) where
  dimension (a,b,c) = dimension a + dimension b + dimension c


type FScalarBasisSpace a = (FinDimVecSpace a, Floating (Scalar a), Field (Scalar a))