module QuadSieve where

-- References
-- https://www.cs.virginia.edu/crab/QFS_Simple.pdf
 

import Math.NumberTheory.Primes (nextPrime, precPrime,Prime,unPrime)
import Math.NumberTheory.Powers.Squares (integerSquareRoot)
import Math.NumberTheory.Moduli.Sqrt (sqrtsModPrime)
import qualified Data.Vector as Vec
import Data.Vector ((!),(//))
import qualified Data.Map.Strict as M
import Data.Maybe (isJust,fromJust)

--- ========================= QUADRATIC SIEVE ============================== --

--- Look for the factors of an odd composite integer n ---
quadSieve :: Integer -> Maybe (Integer,Integer)
quadSieve = undefined

--- ========================== QUADRATIC FUNCTION ========================== --

data Quad = Quad {
  root  :: Integer
, value :: Integer 
} deriving (Eq,Show)

quadratic :: Integer -> Integer -> Quad
quadratic x n = Quad x (x^2 - n)

--- ======================== OPTIMAL B and OPTIMAL M ======================= --

optimalB :: Integer -> Integer
optimalB n = ceiling (optimalB' n)

optimalB' :: Integer -> Double
optimalB' n = exp (p1*p2)
  where 
    n' = fromInteger n
    p1 = sqrt $ log n' * log (log n')
    p2 = 1/2

optimalM :: Integer -> Integer
optimalM n = optimalB n ^ 3

-- ============================= Find Bsmooth Quadratics =============================== --

-- Find b-smooth quadratics x^2 -n  for x in range [0..m]
--- For each b-smooth quadratic return its (x,Q(x) and its prime factorisation)

findBSmoothQuad :: Integer -> Integer -> Integer ->  Vec.Vector (Quad,[(Prime Integer, Integer)])
findBSmoothQuad n b m = 
  let 
    primes = [nextPrime 2 .. precPrime b]
    quads = Vec.fromList [quadratic x n | x <- [0..m]]
    candidates = sieve n primes quads 
    f (q, m) = (q,M.toList m)
  in 
    Vec.map f candidates
   
-- =================================== SIEVE ========================================== --

sieve :: Integer -> [Prime Integer] -> Vec.Vector Quad -> Vec.Vector (Quad, M.Map (Prime Integer) Integer)
sieve n ps quads = 
  let 
    --- Use for dividing ---
    quads' :: Vec.Vector Integer
    quads' = Vec.map value quads
    
    --- Use to keep track of prime factors ---
    emptyFactors :: Vec.Vector (M.Map (Prime Integer) Integer)
    emptyFactors = Vec.replicate (Vec.length quads) M.empty

    (divided,primeFactorization) = Vec.unzip $ sieve' n ps (Vec.zip quads' emptyFactors)

    f :: Integer -> Quad ->  M.Map (Prime Integer) Integer -> Maybe (Quad, M.Map (Prime Integer) Integer)
    f 1 q m = Just (q,m)
    f _ _ _ = Nothing 

    candidates :: Vec.Vector (Maybe (Quad, M.Map (Prime Integer) Integer)) 
    candidates = Vec.zipWith3 f divided quads primeFactorization
  in   
  Vec.map fromJust (Vec.filter isJust candidates)

---- Take a a sequence of primes and sieve ----
sieve' :: Integer -> [Prime Integer] -> Vec.Vector (Integer, M.Map (Prime Integer) Integer) -> Vec.Vector (Integer, M.Map (Prime Integer) Integer)
sieve' n [] vec = vec
sieve' n (p:ps) vec = 
   let 
     s = sieve'' n p  vec
     t = seq s ()
   in 
    sieve' n ps s

---- Take a  primes and sieve ----
sieve'' :: Integer -> Prime Integer -> Vec.Vector (Integer, M.Map (Prime Integer) Integer) -> Vec.Vector (Integer, M.Map (Prime Integer) Integer)
sieve'' n p vec = 
  let 
    --- Find the quadratic residues of x^2 = 0 mod p
    residues :: [Int]
    residues = map fromIntegral (sqrtsModPrime n p)

    p' :: Int
    p' = fromIntegral (unPrime p)

    l :: Int 
    l = Vec.length vec

    --- The set of indices we want to update ---
    indices :: [Int]
    indices = [x| r <- residues, x <- [r, r+p'.. l-1]]   

    --- Values we want to update ----
    vals :: [(Integer, M.Map (Prime Integer) Integer)]
    vals =  getValues vec indices 
    
    --- Given a prime , (quot , prime factorisation) divide quot by p as many times as possible. ---
    --- Then memo the number of times divided 
    f :: Prime Integer -> (Integer, M.Map (Prime Integer) Integer) -> (Integer, M.Map (Prime Integer) Integer)
    f p (n,m) = 
      let 
        p' = unPrime p
        (n', divs) = divideByP n p'
      in 
        (n' , M.insert p divs m)
    
    updatedVals :: [(Integer, M.Map (Prime Integer) Integer)]
    updatedVals =  map (f p) vals

    updates :: [(Int,(Integer, M.Map (Prime Integer) Integer))]
    updates = zip indices updatedVals 

  in 
     vec // updates


--- ================ Helper functions ============================ ---

getValues :: Vec.Vector a -> [Int] -> [a]
getValues _ [] = []
getValues vec (i:is) = (vec ! i) : getValues vec is

divideByP ::  Integer -> Integer -> (Integer, Integer)
divideByP n p = divideByP' n p 0
   
divideByP' :: Integer -> Integer ->Integer -> (Integer, Integer)
divideByP' n p t
  | n == 0         = (n,0)
  | n `mod` p == 0 =  divideByP' (n `div` p) p (t+1)
  | otherwise      = (n,t) 
