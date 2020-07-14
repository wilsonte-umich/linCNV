
# TODO:  modify mark_cells to include Major Variant button
# use it to call a first, obvious cell subset pattern, e.g. cancer

# codify cell lineage sets
cellSetPrototype <- list(
    level = integer, 1==all cells, 2== subset of all cells, 3== subset of level 2 set, ..., changeable
    cells = integer(),
    cnvs  = list()
)

# initialize the top-level set representing all cells
cellSets <- list(
    list( 
        level = 1, # will never change for the top-level set (can't be a subset of anything else)
        cells = cell$Is_accepted,
        cnvs  = list()
    )
)

# if user marked a major variant, add it next during initialization
if(user marked a major variant){
    cellSets <- list(cellSets, list(
        level = 2,
        cells = cell$Is_accepted[mark == majorVariantMark],
        cnvs <- list()
    ))
}

# TODO: when a new cellSet is identified/nominated
newSet <- nominated by Shiny UI, or code? # must set cnvs as 1-item list and cells
newRelationships <- getRelationships(newSet, cellSets)

identicalSetI <- which(newRelationships == NEW_IS_IDENTICAL) # will only ever be one exact match
if(length(identicalSetI) > 0){ # we know this set, just add the new CNV to the list that defines the set
    cellSets[[identicalSetI]]$cnvs <- list(cellSets[[identicalSetI]]$cnvs, newSet$cnvs)
} else { # a new set, place it in the tree
    parentSetIs <- which(newRelationships == NEW_IS_SUBSET) # will always include at least the set of all cells
    parentSetLevels <- sapply(parentSetIs, function(i) cellSets[[i]]$level)
    parentSetI <- which.max(parentSetLevels) # the smallest set that contains newSet

    childSetIs <- which(newRelationships == NEW_IS_SUPERSET) # may be none
    if(length(childSetIs) > 0){
        childSetLevels <- sapply(childSetIs, function(i) cellSets[[i]]$level)
        newSet$level <- min(childSetLevels) # insert new set at top of tree of children and demote the children
        # newSet$level should == cellSets[[parentSetI]]$level + 1
        for(i in childSetIs) cellSets[[i]]$level <- cellSets[[i]]$level + 1 
    }  else {
        newSet$level <- cellSets[[parentSetI]]$level + 1
    }
    cellSets <- list(cellSets, newSet)
}

# function to determine the relationship of a pending cellSet to all others
NEW_IS_IDENTICAL <- 1
NEW_IS_SUBSET    <- 2
NEW_IS_SUPERSET  <- 3
NEW_IN_CONFLICT  <- 4
NEW_IS_UNRELATED <- 5
getRelationships <- Vectorize(function(newSet, oldSet){
    nNewInOld <- sum(newSet$cells %in% oldSet$cells)
    nOldInNew <- sum(oldSet$cells %in% newSet$cells) 
    newIsSubset    <- nNewInOld == length(newSet$cells)
    newIsSuperset  <- nOldInNew == length(oldSet$cells)
    if(newIsSubset & newIsSuperset){
        NEW_IS_IDENTICAL
    } else if(newIsSubset){
        NEW_IS_SUBSET
    } else if(newIsSuperset){
        NEW_IS_SUPERSET    
    } else if(nNewInOld != 0){ # some cells in new are not in old, AND some cells in old are not in new!
        NEW_IN_CONFLICT
    } else {
        NEW_IS_UNRELATED
    }
})

