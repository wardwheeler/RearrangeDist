{- |
Module      :  rearrangeDist.hs 
Description :  Progam to calculate chromosomal rearangement distances
               input fastc file
Copyright   :  (c) 2020 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
License     :  

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.

Maintainer  :  Ward Wheeler <wheeler@amnh.org>
Stability   :  unstable
Portability :  portable (I hope)

-}

module Main where

import System.IO
import System.Environment
import Data.List
import Control.Parallel.Strategies
import Control.Concurrent
import System.IO.Unsafe
import Data.Maybe
import Debug.Trace
import qualified SymMatrix as SM
import Data.CSV
import qualified Data.Map.Strict as Map


-- | getNumThreads gets number of COncurrent  threads
{-# NOINLINE getNumThreads #-}
getNumThreads :: Int
getNumThreads = unsafePerformIO getNumCapabilities

-- | fst3 returns first element of triple
fst3 :: (a,b,c) -> a
fst3 (d,_,_) = d
-- | thd3 returns third element of triple
thd3 :: (a,b,c) -> c
thd3 (_,_,f) = f


-- | parseFastc take raw genome data and recodes in list of pairs of (name, [loci])
parseFastc :: String -> ([(String, [String])], [String])
parseFastc inData =
    let inLines' = lines inData
        inLines = filter (/= []) inLines'
        taxonNames = fmap tail $ filter  ((== '>').head) (words inData)
        taxonDataLines = filter  ((/= '>').head) inLines
        taxonData = fmap words taxonDataLines
        locusNames = nub $ fmap (filter (/= '~')) $ concat taxonData --sort $ fmap (filter (/= '~')) $ 
    in
    -- trace ("There are " ++ (show $ length inLines) ++ " input lines in data file")
    (zip taxonNames taxonData, locusNames)

-- | makePairs takes a list of loci and topology and returns adjacent pairs (orderad)
-- if circular adds a front to back pair
makePairs ::  [String] -> Char -> (String, String) -> [(String, String)]
makePairs inLocusList topology headTail =
    if length inLocusList == 1 then
        if topology == 'c' then [headTail]
        else []
    else
        if (head inLocusList) < (head $ tail inLocusList) then 
            (head inLocusList, head $ tail inLocusList) : makePairs (tail inLocusList) topology headTail
        else (head $ tail inLocusList, head inLocusList) : makePairs (tail inLocusList) topology headTail

-- | getBreakPoint takes optology and rearrange cost and indel cost and returns distance 
getBreakPoint :: (String, [String]) -> (String, [String]) -> Char -> Int -> Int -> Int
getBreakPoint (_, seq1Loci) (_, seq2Loci) topology rearrangeCost inDelCost =
    if null seq1Loci || null seq2Loci then error "Null input chromosome"
    else 
        let in1not2 = seq1Loci \\ seq2Loci
            in2not1 = seq2Loci \\ seq1Loci
            toRemove = concat [in1not2, in2not1]
            indelAdjustment = inDelCost * ((length in1not2) + (length in2not1)) 
            seq1' = seq1Loci \\ toRemove
            seq2' = seq2Loci \\ toRemove
            pairs1 = makePairs seq1' topology (head seq1', last seq1')
            pairs2 = makePairs seq2' topology (head seq2', last seq2')
            distPairs = (pairs1 \\ pairs2) ++ (pairs2 \\ pairs1)
            pairCost = (length distPairs) * rearrangeCost
        in
        pairCost + indelAdjustment

-- | getPairDist takes arguments and calls appropriate distance calculator
getPairDist :: Char -> Char -> Int -> Int -> (String, [String]) -> (String, [String])-> Int
getPairDist method topology rearrangeCost inDelCost seq1 seq2  =
    if method == 'b' then
        getBreakPoint seq1 seq2 topology rearrangeCost inDelCost
    else error ("Method (first letter) " ++ (show method) ++ " not implemented")

-- | getPairs creates a matrix of parwise distances
getPairs :: [(String, [String])] -> [(String, [String])] -> Char -> Char -> Int -> Int -> [[Int]]
getPairs rowSeqs columnSeqs method topology rearrangeCost inDelCost =
    if null rowSeqs then []
    else
        let (chunkSize, _) = quotRem (length columnSeqs) getNumThreads
            firstRowSeq = head rowSeqs
            rowDistList = fmap (getPairDist method topology rearrangeCost inDelCost firstRowSeq) columnSeqs `using` parListChunk chunkSize rdeepseq
        in
        rowDistList : getPairs (tail rowSeqs) columnSeqs method topology rearrangeCost inDelCost 

-- | makeStringList take list of list is retuns nice matrix form 
-- add indel (10 x max value) column and row
makeStringList :: [[Int]]-> Int -> Int -> [String]
makeStringList inListList maxDist numStates =
    if null inListList then
        -- add last row
        let lastRow = (replicate numStates (10 * maxDist)) ++ [0]
            rowString = concat $ intersperse (" ") $ (fmap show lastRow)
        in 
        [rowString]
    else 
        let firstRow = (head inListList) ++ [10 * maxDist] 
            stringList = concat $ intersperse (" ") $ (fmap show firstRow)
        in
        stringList : makeStringList (tail inListList) maxDist numStates

-- | pair2String takes pairs and converts to string
pair2String :: (String, String) -> String
pair2String (x,y) = (x ++ " " ++ y)

-- | checkForZero take a seqeuence index and checks in matrix for a zero value, if so
-- returns teh match indix (in Maybe) else nothing.
checkForZero :: Int -> [[Int]] -> Maybe Int
checkForZero seqIndex distMatrix =
    let firstZero = elemIndex 0 $ take seqIndex (distMatrix !! seqIndex)
    in 
    if firstZero == Nothing then Nothing
    else firstZero

-- | getMinimalStatesAndMatrix takes a matrix and "compresses" removing 0 costs 
-- and collapsing states with 0 distance, assigning them to first taxon with 0 distance 
-- newStateList contains teh updated states and original states for each taxon as states are examined for 
-- uniquness (no zero distances)
getMinimalStatesAndMatrix :: [String] -> [[Int]] -> Int -> Int -> [Int] -> [(String, Int, Int)] -> ([(String, Int)],[[Int]])
getMinimalStatesAndMatrix inSeqs distMatrix seqCounter stateCounter deleteList newStateList =
    if null inSeqs then 
        let tempLowerDiag = SM.fromLists distMatrix
            newLowerDiag = SM.deleteRowsAndColumns tempLowerDiag deleteList
            newMatrix = SM.toFullLists newLowerDiag
            stateList = fmap thd3 newStateList
            taxList = fmap fst3 newStateList
            reducedStateList = zip taxList stateList
        in (reducedStateList, newMatrix) --reverse since prepending pairs
    else 
        let firstSeq = head inSeqs 
            zeroDist = checkForZero seqCounter distMatrix
        in
        if zeroDist == Nothing then -- is a new states to keep 
            getMinimalStatesAndMatrix (tail inSeqs) distMatrix (seqCounter + 1) (stateCounter + 1) deleteList (newStateList ++ [(firstSeq, seqCounter, stateCounter)])
        else 
            let zeroState = fromJust zeroDist -- original number of matched taxon
                matchState = thd3 $ newStateList !! zeroState -- updated state number of taxon with zero distance (could be after other matches)
            in
            --trace (show ) 
            getMinimalStatesAndMatrix (tail inSeqs) distMatrix (seqCounter + 1) stateCounter (seqCounter : deleteList) (newStateList ++ [(firstSeq, seqCounter, matchState)])


-- | makePreefixFreeList takes list of number strings and letter string and makes list of prefix free Strings
makePrefixFreeList :: [Char] -> [Char] -> [String]
makePrefixFreeList numberList letterList =
    if null numberList then []
    else 
        let firstNum = [head numberList]
            letterStringList = fmap (:[]) letterList -- charToString letterList
            newRow = fmap (firstNum ++ ) letterStringList
        in
        newRow ++ (makePrefixFreeList (tail numberList) letterList)



-- |  makeMap creates prefix free symbol list for symbol list
-- limited to 520 (10 * 52) -- needs to be a bit more sophisticated so 
-- no limited to 520
makeMap :: [String] -> Map.Map String String
makeMap symbolList =
    if null symbolList then error "Empty symbol List"
    else if length symbolList > 520 then error "Only making codes for up to 520 states"
    else
        let numSymbols = length symbolList
            letterList = ['a'..'z'] ++ ['A'..'Z'] 
            numberList = ['0'..'9']
            prefixFreeList = makePrefixFreeList numberList letterList
            keyPairList = zip symbolList (take numSymbols prefixFreeList)
        in
        trace ("keys : " ++ (show $ length numberList) ++ " numbers " ++ (show $ length letterList)  
            ++ " lettters " ++ (show $ length prefixFreeList) ++ " " ++ (show $ length keyPairList))
        Map.fromList keyPairList


-- | getNewSymbols takes symbol list and map and returns new symbol as String
getNewSymbols :: Map.Map String String -> [String] -> String
getNewSymbols symbolMap symbols =
    if null symbols then []
    else 
        let firstSymbol = Map.lookup (head symbols) symbolMap
        in
        if firstSymbol == Nothing then error ("Can't find symbol " ++ (head symbols) ++ " in map")  
        else 
            (fromJust firstSymbol) ++ " "  ++  (getNewSymbols symbolMap (tail symbols))

-- | getNewSymbols takes a symbol String returns new symbol as String
getNewSymbol :: Map.Map String String -> String -> String
getNewSymbol symbolMap symbol  =
    if null symbol then []
    else 
        let firstSymbol = Map.lookup symbol symbolMap
        in
        if firstSymbol == Nothing then error ("Can't find symbol " ++ symbol ++ " in map")  
        else 
            (fromJust firstSymbol)


-- | main driver
main :: IO ()
main =
  do
    -- Process arguments
    --  csv file first line taxon/leaf names, subsequent lines are distances--must be symmetrical
    args <- getArgs
    if length args < 6 then error "Need at least 4 arguments: input fastc file, method ('breakpoint' or 'inversion'), 'linear' or 'circular' chomosome, cost (integer) of rearrangement, cost of insertion/delettion, whether all taxa have unique states (full/reduced), and output file stub name "
    else hPutStrLn stderr ("\nReading " ++ show (head args) ++ " input file Output to stub.fastc and stub.tcm ")

    let method = head (args !! 1)
    let topology = head (args !! 2)
    let rearrangeCost = read (args !! 3) :: Int
    let inDelCost = read (args !! 4) :: Int
    let stateForm = args !! 5
    let fileStub = args !! 6

    mapM_ (hPutStrLn stderr) (tail args)

    rawData <- readFile $ head args
    let (inSeqs, locusNames) = parseFastc rawData

    hPutStrLn stderr ("There are " ++ (show $ length inSeqs) ++ " taxa")
    hPutStrLn stderr ("There are " ++ (show $ length locusNames) ++ " loci")

    let distMatrix = getPairs inSeqs inSeqs method topology rearrangeCost inDelCost 
    let maxDist = maximum $ fmap maximum distMatrix
    hPutStrLn stderr ("The maximum distance is " ++ (show maxDist))

    let (newStates, newMatrix) = getMinimalStatesAndMatrix (fmap fst inSeqs) distMatrix 0 0 [] []
    hPutStrLn stderr ("There are " ++ (show $ length newMatrix) ++ " unique states")

    let stateListFull = fmap show [0..((length distMatrix) - 1)]
    let stateListReduced = fmap show [0..((length newMatrix) - 1)]

    let symbolMapFull = makeMap stateListFull
    let symbolMapReduced = makeMap stateListReduced

    let symbolStringFull = getNewSymbols symbolMapFull stateListFull 
    let symbolStringReduced = getNewSymbols symbolMapReduced stateListReduced 

    let symbolListFull = fmap (getNewSymbol symbolMapFull) stateListFull 
    -- let symbolListReduced = fmap (getNewSymbol symbolMapFull) stateListReduced 

    -- add large number row and large number column for indels
    let distString = symbolStringFull : makeStringList distMatrix maxDist (length distMatrix)
    let newDistString = symbolStringReduced : makeStringList newMatrix maxDist (length newMatrix)

    let charDefFull = zip (fmap fst inSeqs) symbolListFull
    let newMappedStates = fmap (getNewSymbol symbolMapReduced) $ fmap show $ fmap snd newStates
    let charDefReduced = zip (fmap fst inSeqs) newMappedStates

    let fastcFile = fileStub ++ ".fastc"
    let tcmFile = fileStub ++ ".tcm"
    let csvOutFile = fileStub ++ ".csv"

    fout <- openFile fastcFile WriteMode
    tout <- openFile tcmFile WriteMode
    cout <- openFile csvOutFile WriteMode

    
    -- CSV output of full matrix
    let nameList = fmap fst inSeqs
    let distMatrixStrings = fmap (fmap show) distMatrix
    let csvString = genCsvFile (nameList : distMatrixStrings)
    hPutStrLn cout csvString
    

    if (head stateForm == 'f') then do 
        {mapM_  (hPutStrLn fout) $ fmap pair2String charDefFull; --need prefix free
        -- Add indel values
        mapM_ (hPutStrLn tout) distString;}
    else do 
        {mapM_  (hPutStrLn fout) $ fmap pair2String charDefReduced; --need prefix free
        -- add indel values
        mapM_ (hPutStrLn tout) newDistString;}

    hClose fout
    hClose cout
    hClose tout

    hPutStrLn stderr "Done"
    