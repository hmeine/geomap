#ifndef CELLPYRAMID_HXX
#define CELLPYRAMID_HXX

namespace vigra {

template<class Segmentation, class CellStatistics>
class CellPyramid
{
    struct Level
    {
        unsigned int level;
        Segmentation segmentation;
    };

    typedef typename Segmentation::CellInfo CellInfo;
    typedef typename Segmentation::NodeInfo NodeInfo;
    typedef typename Segmentation::EdgeInfo EdgeInfo;
    typedef typename Segmentation::FaceInfo FaceInfo;
    typedef typename Segmentation::DartTraverser DartTraverser;

    enum OperationType { RemoveIsolatedNode,
                         MergeFaces,
                         RemoveBridge,
                         MergeEdges,
                         // composed Operations:
                         RemoveEdge,
                         RemoveEdgeWithEnds };

    struct Operation
    {
        OperationType  type_;
        DartTraverser  param_;

        Operation(OperationType type, DartTraverser param)
        : type_(type), param_(param) {}
    };

    typedef std::map<unsigned int, std::pair<Segmentation, CellStatistics> >
        CheckpointMap;
    CheckpointMap checkpoints_;

    typedef std::vector<Operation>
        History;
    History history_;

    unsigned int   currentLevel_;
    Segmentation   currentSegmentation_;
    CellStatistics currentCellStatistics_;

    FaceInfo &removeIsolatedNodeInternal(const DartTraverser & dart)
    {
        currentCellStatistics_.preRemoveIsolatedNode(dart);
        FaceInfo &result(currentSegmentation_.removeIsolatedNode(dart));
        currentCellStatistics_.postRemoveIsolatedNode(dart, result);
        return result;
    }

    FaceInfo &mergeFacesInternal(const DartTraverser & dart)
    {
        currentCellStatistics_.preMergeFaces(dart);
        FaceInfo &result(currentSegmentation_.mergeFaces(dart));
        currentCellStatistics_.postMergeFaces(dart, result);
        return result;
    }

    FaceInfo &removeBridgeInternal(const DartTraverser & dart)
    {
        currentCellStatistics_.preRemoveBridge(dart);
        FaceInfo &result(currentSegmentation_.removeBridge(dart));
        currentCellStatistics_.postRemoveBridge(dart, result);
        return result;
    }

    EdgeInfo &mergeEdgesInternal(const DartTraverser & dart)
    {
        currentCellStatistics_.preMergeEdges(dart);
        EdgeInfo &result(currentSegmentation_.mergeEdges(dart));
        currentCellStatistics_.postMergeEdges(dart, result);
        return result;
    }

    CellInfo &performOperation(Operation &op)
    {
        switch(op.type_)
        {
          case RemoveIsolatedNode:
          {
              return removeIsolatedNodeInternal(op.param_);
          }
          case MergeFaces:
          {
              return mergeFacesInternal(op.param_);
          }
          case RemoveBridge:
          {
              return removeBridgeInternal(op.param_);
          }
          case MergeEdges:
          {
              return mergeEdgesInternal(op.param_);
          }
          case RemoveEdge:
          {
              return (op.param_.leftFaceLabel() == op.param_.rightFaceLabel() ?
                      removeBridgeInternal(op.param_) :
                      mergeFacesInternal(op.param_));
          }
          case RemoveEdgeWithEnds:
          {
              EdgeInfo &removedEdge = currentSegmentation_.edge(op.param_.edgeLabel());
              bool removedEdgeIsLoop = removedEdge.start.startNodeLabel() ==
                                       removedEdge.end.startNodeLabel();

              FaceInfo &result = (op.param_.leftFaceLabel() == op.param_.rightFaceLabel() ?
                                  removeBridgeInternal(op.param_) :
                                  mergeFacesInternal(op.param_));

              if(removedEdge.start.recheckSingularity())
                  removeIsolatedNodeInternal(removedEdge.start);
              else
              {
                  // recheckSingularity() modified removedEdge.start
                  // to point to a valid edge pixel again, so it is
                  // not really removedEdge.start anymore.. ;-)
                  DartTraverser changedNode(removedEdge.start);
                  // test if the DartTraverser has degree 2 and the
                  // edge is no loop:
                  if((changedNode.nextSigma() != removedEdge.start) &&
                     (changedNode.edgeLabel() != removedEdge.start.edgeLabel()) &&
                     (changedNode.nextSigma() == removedEdge.start))
                  {
                      mergeEdgesInternal(changedNode);
                  }
              }

              if(!removedEdgeIsLoop)
              {
                  if(removedEdge.end.recheckSingularity())
                      removeIsolatedNodeInternal(removedEdge.end);
                  else
                  {
                      DartTraverser changedNode(removedEdge.end);
                      if((changedNode.nextSigma() != removedEdge.end) &&
                         (changedNode.edgeLabel() != removedEdge.end.edgeLabel()) &&
                         (changedNode.nextSigma() == removedEdge.end))
                      {
                          mergeEdgesInternal(changedNode);
                      }
                  }
              }

              return result;
          }
        }

        vigra_fail("Unknown operation type in CellPyramid<>::performOperation!");
        return currentSegmentation_.face(0);
    }

    CellInfo &gotoLevel(unsigned int level)
    {
        vigra_precondition(level < levelCount(),
                           "trying to access non-existing level (too high index)!");

        if(currentLevel_ > level)
        {
            typename CheckpointMap::iterator lastCheckpointIt =
                checkpoints_.upper_bound(level);
            --lastCheckpointIt;
            currentLevel_ = lastCheckpointIt->first;
            currentSegmentation_ = lastCheckpointIt->second.first;
            currentCellStatistics_ = lastCheckpointIt->second.second;
        }

        while(currentLevel_ < level)
        {
            CellInfo &result= performOperation(history_[currentLevel_]);
            ++currentLevel_;
            if(currentLevel_ == level)
                return result;
        }

        return currentSegmentation_.face(0); // something HAS to be returned... ;-)
    }

  public:
    void storeCheckpoint()
    {
        checkpoints_.insert(std::make_pair(levelCount()-1, std::make_pair(currentSegmentation_, currentCellStatistics_)));
    }

    CellPyramid(const Segmentation &level0,
                const CellStatistics &level0Stats = CellStatistics())
    : currentSegmentation_(level0),
      currentCellStatistics_(level0Stats)
    {
        storeCheckpoint();
    }

    FaceInfo &removeIsolatedNode(const DartTraverser & dart)
    {
        history_.push_back(Operation(RemoveIsolatedNode, dart));
        return static_cast<FaceInfo &>(gotoLevel(levelCount()-1));
    }

    FaceInfo &mergeFaces(const DartTraverser & dart)
    {
        history_.push_back(Operation(MergeFaces, dart));
        return static_cast<FaceInfo &>(gotoLevel(levelCount()-1));
    }

    FaceInfo &removeBridge(const DartTraverser & dart)
    {
        history_.push_back(Operation(RemoveBridge, dart));
        return static_cast<FaceInfo &>(gotoLevel(levelCount()-1));
    }

    EdgeInfo &mergeEdges(const DartTraverser & dart)
    {
        history_.push_back(Operation(MergeEdges, dart));
        return static_cast<EdgeInfo &>(gotoLevel(levelCount()-1));
    }

    FaceInfo &removeEdge(const DartTraverser & dart)
    {
        history_.push_back(Operation(RemoveEdge, dart));
        return static_cast<FaceInfo &>(gotoLevel(levelCount()-1));
    }

    FaceInfo &removeEdgeWithEnds(const DartTraverser & dart)
    {
        history_.push_back(Operation(RemoveEdgeWithEnds, dart));
        return static_cast<FaceInfo &>(gotoLevel(levelCount()-1));
    }

    const Segmentation &currentLevel() const
    {
        return currentLevel_;
    }

    const Segmentation &currentSegmentation() const
    {
        return currentSegmentation_;
    }

    const CellStatistics &currentCellStatistics() const
    {
        return currentCellStatistics_;
    }

    const Segmentation &segmentation(unsigned int level)
    {
        gotoLevel(level);
        return currentSegmentation_;
    }

    void cutHead()
    {
        history_.erase(history.begin() + currentLevel_, history.end());
        checkpoints_.erase(checkpoints_.upper_bound(currentLevel_),
                           checkpoints_.end());
    }

    unsigned int levelCount() const
    {
        return history_.size() + 1;
    }
};

} // namespace vigra

#endif // CELLPYRAMID_HXX
