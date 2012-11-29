{-# LANGUAGE ConstraintKinds        #-}
{-# LANGUAGE FlexibleContexts       #-}
{-# LANGUAGE TypeFamilies           #-}
{-# LANGUAGE FlexibleInstances      #-}
{-# LANGUAGE ScopedTypeVariables    #-}
{-# LANGUAGE UndecidableInstances   #-}


module Data.Uncertain( FScalarBasisSpace
                     , Approximate(..)
                     , exactly
                     , (+/-), (±), (+-)
                     )  where


import Data.Uncertain.VectorSpace

import Data.Packed.Matrix
import Numeric.LinearAlgebra


import Control.Functor.Constrained 
import Control.Applicative.Constrained

import Control.Applicative

import Data.List





-- | 'Approximate' values represent measurements that have either been taken
-- on some physical system, or are result of some kind of approximate computation.
-- On value may (and in general will) consist of multiple quantities, stored
-- in a 'VectorSpace' container.
data Approximate a = Approximate
       { expectValue :: a      -- ^ estimated / central value of the error distribution. @a@ is some /n/-dimensional vector space.
       , uncertainSpan :: [a]  -- ^ at most /n/ ortho/gon/al vectors that, with std-normal-distributed coefficients, represent the error distribution.
       }


infixl 6 +/-   -- like '+'
infixl 6 +-

-- | The equivalent '(+/-)' and '(±)' do the obvious thing: combining an \"exact\" value
-- with an uncertainty annotation of the same type, like in @23 +/- 2@.
(+/-) :: FScalarBasisSpace a => a -> a -> Approximate a
v +/- u = Approximate v [u]

(±) :: FScalarBasisSpace a => a -> a -> Approximate a
v ± u = Approximate v [u]

-- | '(+-)' adds more uncertainty-base spanning vectors. For instance, a measurement
-- result consisting of a position /x/ with uncertainty /σx/ and a velocity /vₓ ± σvₓ/
-- might be written @(x,v) ± (σx,0)+-(0,σv)@, which looks admittedly cumbersome compared
-- to simply @(x ± σx, v ± σv)@. However, this way of writing is more general: it can
-- account for arbitrary /correlations/ between the quantities' measurement errors.
(+-) :: FScalarBasisSpace a => Approximate a -> a -> Approximate a
Approximate v us +- u = Approximate v (u:us)


exactly :: FScalarBasisSpace a => a -> Approximate a
exactly v = Approximate v []



-- | The 'Num' and 'Fractional' instances of 'Approximate' only use standard
-- (\"Gaussian\") error propagation formulas, which is equivalent (but more
-- CPU-efficient) to the general PCA-based propagation (e.g. @a+b@ vs.
-- @cfmap(+) a b@ in the case of addition and subtraction, and also becomes
-- equivalent for the other primitive operations when the uncertainties are
-- sufficiently small.
-- These instances can, however, not be used for more complex Functions in
-- which variables might get used more than once, because the errors in
-- seperate `Approximate' values are always assumed as /independent/, which
-- two uses of the same measured/calculated value obviously aren't.

instance (RealFloat a, FScalarBasisSpace a) => Num (Approximate a) where
  
  fromInteger = exactly . fromInteger
  
  (Approximate a []) + (Approximate b us) = Approximate (a+b) us
  (Approximate a us) + (Approximate b []) = Approximate (a+b) us
  (Approximate a [ua]) + (Approximate b [ub])
      = a+b +/- (sqrt $ ua^2 + ub^2)
     
  (Approximate a []) - (Approximate b us) = Approximate (a-b) us
  (Approximate a us) - (Approximate b []) = Approximate (a-b) us
  (Approximate a [ua]) - (Approximate b [ub])
      = a-b +/- (sqrt $ ua^2 + ub^2)
  
  (Approximate a []) * (Approximate b us) = Approximate (a*b) $ map (a*) us
  (Approximate a us) * (Approximate b []) = Approximate (a*b) $ map (b*) us
  (Approximate a [ua]) * (Approximate b [ub])
      = a*b +/- sqrt((a*ub)^2 + (b*ua)^2)
  
  negate (Approximate a us) = Approximate (-a) us
  
  abs (Approximate a us) = Approximate (abs a) us
  
  signum (Approximate a us) = exactly $ signum a
  

instance (RealFloat a, FScalarBasisSpace a) => Fractional (Approximate a) where
  
  fromRational = exactly . fromRational
  
  recip (Approximate a []) = exactly $ recip a
  recip (Approximate a [u]) = recip a +/- u / a^2
  
  (Approximate a []  ) / (Approximate b []  ) = exactly $ a/b
  (Approximate a [u] ) / (Approximate b []  ) = a/b +/- u/b
  (Approximate a []  ) / (Approximate b [u] ) = a/b +/- a*u/(b^2)
  (Approximate a [ua]) / (Approximate b [ub])
      = a/b +/- sqrt((a*ub/b^2)^2 + (ua/b)^2)
  


instance (Show a) => Show (Approximate a) where
  showsPrec _ (Approximate v []) = ("exactly "++) . showsPrec 9 v
  showsPrec fxty x@(Approximate v us)
    | fxty>=6    = ('(':) . showsPrec 0 x . (')':)
    | otherwise  = showsPrec 6 v . (" +/- "++) . foldr1 (.)
                      ( intersperse (" +- "++) $ map (showsPrec 6) us )




-- | General uncertainty propagation, combining propagated "unit errors", by means 
-- of a left singular value decomposition, to a resultant set of basis errors.
-- It's essentially principal component analysis.
-- 
-- Disclaimer: this method has not yet been thoroughly tested; it does work for
-- simple examples but should still be tested for complicated multi-dimensional
-- error propagation problems by comparing to a monte-carlo simulation.
reduceEllipsoidRelevantSpan :: forall a. FScalarBasisSpace a => [a] -> [a]
reduceEllipsoidRelevantSpan [] = []
reduceEllipsoidRelevantSpan (v:vs) = vs'
 where vDecomp :: [Scalar a]  
       (basis,vDecomp) = unzip $ decompose v
       decomps :: [[Scalar a]]
       decomps = vDecomp : map( (`map`basis) . decompose' ) vs
       lSVD :: Matrix (Scalar a)
       (lSVD, svs) = leftSV . (length basis >< inDim)
                       . concat $ transpose decomps
       pqd_vs' :: [Vector (Scalar a)]
       pqd_vs' = map (uncurry $ scale . sclCorrect) . filter ((>0) . fst)
                   $ zip (toList svs) (toColumns lSVD)
       vs' = map (recompose . zip basis . toList) pqd_vs'
       
       sclCorrect = realToFrac
       
       inDim = length vs + 1


-- | 'Approximate' is an instance of 'CFunctor' as well as 'CApplicative',
-- which allows any analytical computation to be performed on uncertain
-- values and automatic uncertainty propagation to be performed.
instance CFunctor Approximate where
  type CFunctorCtxt Approximate a = FScalarBasisSpace a
  f `cfmap` Approximate v us = Approximate v' us'
   where v' = f v
         us' = map(^/sqrt(2)) . reduceEllipsoidRelevantSpan
              . map ((v' ^-^) . f) $ [(v ^-^), (v ^+^)] <*> us

instance CApplicative Approximate where
  type CApplicativeCtxt Approximate a = ()
  cpure = exactly
  cliftA2 ifx (Approximate a uas) (Approximate b ubs) = Approximate v uvs
   where v = a `ifx` b
         uvs = map(^/sqrt(2)) . reduceEllipsoidRelevantSpan 
              $ map (v ^-^) =<<
                      [ map (`ifx` b)
                           ([(a ^-^), (a ^+^)] <*> uas)
                      , map (a `ifx`)
                           ([(b ^-^), (b ^+^)] <*> ubs) ]



 