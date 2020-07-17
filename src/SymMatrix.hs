{- |
Module      :  SymMatrix.hs 
Description :  Progam to manipulate square symmetric lower diagonal matrices with diagnonal values
                as if they were normal matrices
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

module SymMatrix (empty, dim, fromLists, Matrix, 
                   SymMatrix.null, cols, rows,
                   (!), toLists, toRows, fromRows,
                   toFullLists, getFullRow,
                   isSymmetric, updateMatrix, 
                   unsafeUpdateMatrix,
                   addMatrixRow, addMatrices,
                   deleteRowsAndColumns, showMatrixNicely) where

-- import Debug.Trace
import qualified Data.List as L
import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import qualified Data.Sort as S

-- | Matrix type as Vector of Vectors
type Matrix a = V.Vector (V.Vector a)

-- | functions for triples
fst3 :: (a,b,c) -> a
fst3 (d,_,_) = d

-- | empty matrix value from Vector
empty :: Matrix a
empty = G.empty

-- | dim returns dimmension (rows = cols)
dim :: (Eq a) => Matrix a -> (Int, Int)
dim inMatrix = 
    if SymMatrix.null inMatrix then (0,0) 
    else (V.length inMatrix, V.length inMatrix)

-- | rows returns number rows in Matrix (rows = cols)
rows :: (Eq a) => Matrix a -> Int
rows inM = fst $ dim inM

-- | cols returns number cols in Matrix (rows = cols)
cols :: (Eq a) => Matrix a -> Int
cols inM = fst $ dim inM

-- | null returns True of row number is 0
null :: (Eq a) => Matrix a -> Bool
null inMatrix = inMatrix == empty

-- | isSymmetric is true by defineition--when creted error if not
isSymmetric :: (Eq a) => Matrix a -> Bool
isSymmetric inM = 
    if SymMatrix.null inM then error "Nulll martix in isSymmetric"
    else True

-- | fromRows creates a lower diagonal matrix (with diagonal)
fromRows :: (Eq a, Show a) => [V.Vector a] -> Matrix a
fromRows inVectList = fromLists $ fmap V.toList inVectList

-- | toRows converts a Matrix to a list of Vectors
-- unequal in length
toRows :: (Eq a, Show a) => Matrix a -> [V.Vector a]
toRows inM = V.toList inM 

-- | fromLists takes list of list of a and returns lower diagnoal (with diagonal)
-- matrix as Matrix (Vector of Vectors)
fromLists :: (Eq a, Show a) => [[a]] -> Matrix a
fromLists inListList =
    if L.null inListList then empty
    else 
        let initialSquare = V.map V.fromList $ V.fromList inListList
            colsH = V.length $ V.head initialSquare
            rowsH = V.length initialSquare
        in
        if (colsH /= rowsH) then error ("Input matrix is not square " ++ show inListList)
        else 
            let indexPairs = cartProd [0..(rowsH - 1)] [0..(rowsH - 1)]
                sym = checkSymmetry initialSquare indexPairs
            in
            if not sym then error "Input matrix not symmetrical"
            else makeLowerDiag initialSquare 0 rowsH

-- | toLists takes a Matrix and returns a list of lists (not all same length)
toLists :: (Eq a, Show a) => Matrix a -> [[a]]
toLists inM =
    if SymMatrix.null inM then []
    else 
        V.toList $ V.map V.toList inM

-- | toFullLists takes a Matrix and returns a list of lists of full length
-- square matrix
toFullLists :: (Eq a, Show a) => Matrix a -> [[a]]
toFullLists inM =
    if SymMatrix.null inM then []
    else 
        fmap (getFullRow inM) [0..((rows inM) - 1)]

-- | getFullRow returns a specific full row (is if matrix were square)
-- as a list
getFullRow :: (Eq a, Show a) => Matrix a -> Int -> [a]
getFullRow inM index =
    if SymMatrix.null inM then []
    else
        let firstPart = V.toList $ inM V.! index -- initial [0..index] of elements
            restMatrix = V.drop (index + 1) inM
            restByColumn = V.toList $ V.map (V.! (index + 1)) restMatrix
        in
        firstPart ++ restByColumn


-- | indexing lower diag matrix
(!) :: Matrix a -> (Int, Int) -> a 
(!) inM (iIndex,jIndex) =
    if iIndex > jIndex then 
        (inM V.! iIndex) V.! jIndex
    else
        (inM V.! jIndex) V.! iIndex

-- | makeLowerDiag take a Vector of Vetors (Matrix) and retusn a lower diagonal matrix
-- including diagonal
makeLowerDiag :: (Eq a) => Matrix a -> Int -> Int -> Matrix a
makeLowerDiag inM row numRows=
    if SymMatrix.null inM then error "Input matrix is empty in makeLowerDiag"
    else if (row == numRows) then V.empty
    else 
        let origRow = inM V.! row
            newRow = V.take (row + 1) origRow
        in
        V.cons newRow  (makeLowerDiag inM (row + 1) numRows)

-- | cartesian product of two lists 
cartProd :: [a] -> [a] -> [(a,a)]
cartProd xs ys = [(x,y) | x <- xs, y <- ys]

-- | checkSymmetry takes a Vector of VEctors of a 
-- and checks for symmetry
checkSymmetry :: (Eq a, Show a) => Matrix a -> [(Int, Int)] -> Bool
checkSymmetry inVV pairList =
    if SymMatrix.null inVV then error ("Input matrix is empty")
    else if L.null pairList then True
    else
        let (iIndex, jIndex) = head pairList
            firstCheck = ((inVV V.! iIndex) V.! jIndex) == ((inVV V.! jIndex) V.! iIndex)
        in
        if firstCheck then checkSymmetry inVV (tail pairList)
        else error ("Matrix is not symmetrical:" ++ (show (iIndex, jIndex)) ++ "=>" ++ (show ((inVV V.! iIndex) V.! jIndex)) ++ " /= " ++ (show ((inVV V.! jIndex) V.! iIndex)))

-- | addMatrixRow add a row to existing matrix to extend Matrix dimension
-- used when adding HTU distances to existing distance matrix as Wagner tree is built
addMatrixRow :: (Eq a) => Matrix a -> V.Vector a -> Matrix a
addMatrixRow inM newRow =
    if V.null newRow then inM
    else 
        inM `V.snoc` newRow


-- | addMatrices  adds a Matrix to existing matrix to extend Matrix dimension
addMatrices :: (Eq a) => Matrix a -> Matrix a -> Matrix a
addMatrices inM newMatrix =
    if SymMatrix.null newMatrix then inM
    else 
        inM V.++ newMatrix

-- | reIndexTriple taske (i,j,k) and returns (max i j, min i j, k)
reIndexTriple :: (Ord a) => (a, a, b) -> (a, a, b)
reIndexTriple trip@(iIndex, jIndex, value) =
    if iIndex > jIndex then trip
    else (jIndex, iIndex, value)

-- | updateMatrix takes a list of triples and update matrix
-- update all at once checking for bounds
-- could naively do each triple in turn, but would be alot of copying
updateMatrix :: (Eq a, Show a, Ord a) => Matrix a -> [(Int, Int, a)] -> Matrix a
updateMatrix inM modList =
    if L.null modList then inM
    else 
        let orderedTripleList = S.uniqueSort $ fmap reIndexTriple modList
            minRow = fst3 $ head orderedTripleList
            maxRow = fst3 $ last orderedTripleList
        in
        if minRow < 0 then error ("Update matrix out of bounds: " ++ (show orderedTripleList))
        else if maxRow >= rows inM then error ("Update matrix out of bounds, row = " ++ (show $ rows inM) ++ " and trying to update row " ++ (show maxRow))
        else
            let firstPart = V.unsafeTake minRow inM
                restPart  = V.unsafeDrop minRow inM
                modifiedRemainder = updateRows restPart orderedTripleList minRow 
            in
            -- trace ("Modifying : " ++ (show modList) ++ "\n" ++ showMatrixNicely firstPart) (
            -- trace (showMatrixNicely $ addMatrices firstPart modifiedRemainder) 
            addMatrices firstPart modifiedRemainder
            -- )

-- | unsafeUpdateMatrix unsafe version of updateMatrix
unsafeUpdateMatrix :: (Eq a, Show a, Ord a) => Matrix a -> [(Int, Int, a)] -> Matrix a
unsafeUpdateMatrix inM modList =
    if L.null modList then inM
    else 
        let orderedTripleList = S.uniqueSort $ fmap reIndexTriple modList
            minRow = fst3 $ head orderedTripleList
            firstPart = V.unsafeTake minRow inM
            restPart  = V.unsafeDrop minRow inM
            modifiedRemainder = updateRows restPart orderedTripleList minRow 
        in
        addMatrices firstPart modifiedRemainder
        
-- | updateRows takes the section of the matrix containing rows that wil be modified
-- (some not) and modifes or copies rows and rerns a Matrix (vector of roow vectors)
updateRows :: (Show a, Eq a) => Matrix a -> [(Int, Int, a)] -> Int -> Matrix a
updateRows inM tripList currentRow =
    if L.null tripList then inM
    else if SymMatrix.null inM then error ("Matrix is empty and there are modifications that remain: " ++ (show $ tripList))
    else 
        let (rowIndex, columnIndex, value) = L.head tripList
            firstOrigRow = V.head inM
        in
        if currentRow /= rowIndex then firstOrigRow `V.cons` (updateRows (V.tail inM) tripList (currentRow + 1))
        else -- account for multiple modifications to same row
            let (newRow, newTripList) = modifyRow firstOrigRow columnIndex value currentRow (L.tail tripList)
            in
            -- This for debug--remove after test
            --if (V.length newRow) /= (V.length firstOrigRow) then error ("Modified row not correct length " ++ (show newRow) ++ " -> " ++ (show firstOrigRow))
            --else 
            newRow `V.cons` (updateRows (V.tail inM) newTripList (currentRow + 1))


-- | modifyRow takes an initial modification (column and value) and then checks to see if there are more modifications in that
-- row (rowNumber) in the remainder of the list of modifications, returning the new row and mod list as a pair
-- assumes that sorted triples sort by first, second, then third elements
modifyRow :: V.Vector a -> Int -> a -> Int -> [(Int, Int, a)] -> (V.Vector a, [(Int, Int, a)])
modifyRow inRow colIndex value rowNumber modList =
    if colIndex >= (V.length inRow) then error ("Column to modify is outside length of row " ++ (show $ (rowNumber, colIndex)))
    else 
        let firstPart = V.unsafeTake colIndex inRow
            remainderPart = V.unsafeDrop (colIndex + 1) inRow
            newRow = firstPart V.++ (value `V.cons` remainderPart)
        in 
        if L.null modList then (newRow, modList)
        else continueRow (firstPart `V.snoc` value) inRow (colIndex + 1) rowNumber modList
            

-- | continueRow continues to modify a row with multiple column modifcations
-- assumes that sorted triples sorted by first, second, then third elements
continueRow :: V.Vector a ->V.Vector a -> Int -> Int -> [(Int, Int, a)] -> (V.Vector a, [(Int, Int, a)])
continueRow partRow origRow colIndex rowNumber modList =
    if colIndex == V.length origRow then (partRow, modList) --completed row
    else if L.null modList then                             --no more modifications
        (partRow V.++ (V.unsafeDrop colIndex origRow), modList)
    else 
        let (nextRowNumber, nextColIndex, nextValue) = L.head modList
        in
        if nextRowNumber /= rowNumber then (partRow V.++ (V.unsafeDrop colIndex origRow), modList)
        else 
            if nextColIndex /= colIndex then continueRow (partRow `V.snoc` (origRow V.! colIndex)) origRow (colIndex + 1) rowNumber modList
            else continueRow (partRow `V.snoc` nextValue) origRow (colIndex + 1) rowNumber (L.tail modList)

-- | makeNiceRow pretty preints a list
makeNiceRow :: (Show a) => V.Vector a -> String
makeNiceRow aVect =
  if V.null aVect then "\n"
  else
    show (V.head aVect) ++ " " ++  makeNiceRow (V.tail aVect)

-- | showNicely pretty prins matrix
showMatrixNicely :: (Show a, Eq a) => Matrix a -> String
showMatrixNicely inM =
  let mRows = rows inM
      mCols = cols inM
      niceRows = V.map makeNiceRow inM
  in
  ("Dimensions: :" ++ show mRows ++ " " ++ show mCols ++ "\n" ++ concat niceRows)

-- | deleteRowsAndColumns take a list of rows (and same index for columns)
-- to delete from Matrix. Uses lisyt to do in single pass
deleteRowsAndColumns :: (Show a, Eq a) => Matrix a -> [Int] -> Matrix a
deleteRowsAndColumns inM deleteList =
    if L.null deleteList then inM
    else deleteRC inM deleteList (rows inM) 0


-- | deleteRC takes matri delete list and counter to delte coumns and rows
deleteRC :: (Show a, Eq a) => Matrix a -> [Int] -> Int -> Int -> Matrix a
deleteRC inM deleteList origRows rowCounter =
    if rowCounter == origRows then empty
    else 
        let firstRow = V.head inM
            toKeep = rowCounter `L.notElem` deleteList
            newRow = deleteColumn firstRow deleteList (rowCounter + 1) 0
        in
        if toKeep then newRow `V.cons` (deleteRC (V.tail inM) deleteList origRows (rowCounter + 1))
        else deleteRC (V.tail inM) deleteList origRows (rowCounter + 1)

-- | deleteColumn takes a row of a matrix (lower diagnonal), its length, 
-- a list of cilumns to delete and a column counter and creates a new row
deleteColumn :: (Show a, Eq a) => V.Vector a -> [Int] -> Int -> Int -> V.Vector a
deleteColumn origRow deleteList rowLength colCounter =
    if colCounter == rowLength then V.empty
    else 
        let firstValue = V.head origRow 
            toKeep = colCounter `L.notElem` deleteList
        in
        if toKeep == True then firstValue `V.cons` (deleteColumn (V.tail origRow) deleteList rowLength (colCounter + 1))
        else deleteColumn (V.tail origRow) deleteList rowLength (colCounter + 1)