{-# LANGUAGE ConstraintKinds        #-}
{-# LANGUAGE FlexibleContexts       #-}
{-# LANGUAGE TypeFamilies           #-}
{-# LANGUAGE FlexibleInstances      #-}
{-# LANGUAGE ScopedTypeVariables    #-}
{-# LANGUAGE UndecidableInstances   #-}
{-# LANGUAGE GADTs                  #-}


module Data.Uncertain( where


import Data.Uncertain.VectorSpace


import Control.Functor.Constrained 
import Control.Applicative.Constrained

import Control.Applicative

import Data.List




-- | General uncertainty propagation, combining propagated "unit errors", by means 
-- of a left singular value decomposition, to a resultant set of basis errors.
-- It's essentially principal component analysis.
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
       
       sclCorrect = realToFrac -- . (*sqrt(fromIntegral outDim/fromIntegral inDim) )
       
       (inDim, outDim) = (length vs + 1, dimension v)


infixl 6 +/-
infixl 6 +-
data Uncertain a = Uncertain
       { expectValue :: a      -- ^ estimated / central value of the error distribution. @a@ is some /n/-dimensional vector space.
       , uncertainSpan :: [a]  -- ^ at most /n/ ortho/gon/al vectors that, with std-normal-distributed coefficients, represent the error distribution.
       }


(+/-) :: FScalarBasisSpace a => a -> a -> Approximate a
v +/- u = Approximate v [u]
(+-) :: FScalarBasisSpace a => Approximate a -> a -> Approximate a
Approximate v us +- u = Approximate v (u:us)

exactly :: FScalarBasisSpace a => a -> Approximate a
exactly v = Approximate v []



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

instance CFunctor Approximate where
  type CFunctorCtxt Approximate a = FScalarBasisSpace a
  f `cfmap` Approximate v us = Approximate v' us'
   where v' = f v
         us' = reduceEllipsoidRelevantSpan
              . map ((v' ^-^) . f) $ [(v ^-^), (v ^+^)] <*> us

instance CApplicative Approximate where
  type CApplicativeCtxt Approximate a = ()
  cpure = exactly
  cliftA2 ifx (Approximate a uas) (Approximate b ubs) = Approximate v uvs
   where v = a `ifx` b
         uvs = map(^/sqrt(2)) $ reduceEllipsoidRelevantSpan devs
         devs = map (v ^-^) =<<
                      [ map (`ifx` b)
                           ([(a ^-^), (a ^+^)] <*> uas)
                      , map (a `ifx`)
                           ([(b ^-^), (b ^+^)] <*> ubs) ]



 

-- data SmoothFunApproxResult a b = SmoothFunApproxResult
--        { smFunResEstimate :: Approximate b
--        , smFunContinuationDomain :: Approximate a
--        , smFunContinuation :: SmoothFunApprox a b }
-- 
-- newtype SmoothFunApprox a b = SmoothFunApprox
--        { runSmFunApprox :: a -> SmoothFunApproxResult a b }
-- 
-- instance Category SmoothFunApprox where
--   id = SmoothFunApprox (\a -> SmoothFunApproxResult
--                                 ( Approximate a Nothing )
--                                 ( Approximate a Nothing )
--                                 ( id )                    )
--   SmoothFunApprox f . SmoothFunApprox g
--      = SmoothFunApprox (\a -> SmoothFunApproxResult
--                                 ( 