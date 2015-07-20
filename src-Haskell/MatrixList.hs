module MatrixList
( MatrixList
) where

type MatrixList = ( Int, Int, [ (Int, Int, Double) ] )

-- matrix height, width, and list of entries (row, column, value) ;
-- including zero-value entries is optional
