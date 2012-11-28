{-# LANGUAGE TypeFamilies                  #-}
{-# LANGUAGE ConstraintKinds               #-}
{-# LANGUAGE FlexibleInstances             #-}
{-# LANGUAGE UndecidableInstances          #-}


module Data.Uncertain.VectorSpace( module Data.VectorSpace
                                 , module Data.Basis
                                 , FScalarBasisSpace
                                 )  where


import Data.VectorSpace        -- from Vectorspace package
import Data.Basis

import Data.Packed.Matrix      -- from hmatrix package
import Numeric.LinearAlgebra





type FScalarBasisSpace a = (HasBasis a, Floating (Scalar a), Field (Scalar a))