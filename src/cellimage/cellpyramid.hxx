#ifndef CELLPYRAMID_HXX
#define CELLPYRAMID_HXX

#include <vector>
#include <map>

namespace vigra {

template<class SEGMENTATION, class CELLSTATISTICS>
class CellPyramid
{
  public:
    typedef SEGMENTATION Segmentation;
    typedef CELLSTATISTICS CellStatistics;

    typedef typename Segmentation::CellInfo CellInfo;
    typedef typename Segmentation::NodeInfo NodeInfo;
    typedef typename Segmentation::EdgeInfo EdgeInfo;
    typedef typename Segmentation::FaceInfo FaceInfo;
    typedef typename Segmentation::DartTraverser DartTraverser;

    struct Level
    {
        unsigned int   level;
        Segmentation   segmentation;
        CellStatistics cellStatistics;

        Level(unsigned int l,
              const Segmentation &s,
              const CellStatistics &c)
        : level(l), segmentation(s), cellStatistics(c)
        {}

        
    };

  private:
    enum OperationType { RemoveIsolatedNode,
                         MergeFaces,
                         RemoveBridge,
                         MergeEdges,
                         // composed Operations:
                         RemoveEdge,
                         RemoveEdgeWithEnds };

    struct Operation
    {
        OperationType                      type;
        typename DartTraverser::Serialized param;

        Operation(OperationType opType, DartTraverser paramDart)
        : type(opType), param(paramDart.serialize()) {}
    };

    typedef std::map<unsigned int, Level>
        CheckpointMap;
    CheckpointMap checkpoints_;

    typedef std::vector<Operation>
        History;
    History history_;

    Level currentLevel_;

    FaceInfo &removeIsolatedNodeInternal(const DartTraverser & dart)
    {
        currentLevel_.cellStatistics.preRemoveIsolatedNode(dart);
        FaceInfo &result(currentLevel_.segmentation.removeIsolatedNode(dart));
        currentLevel_.cellStatistics.postRemoveIsolatedNode(dart, result);
        return result;
    }

    FaceInfo &mergeFacesInternal(const DartTraverser & dart)
    {
        currentLevel_.cellStatistics.preMergeFaces(dart);
        FaceInfo &result(currentLevel_.segmentation.mergeFaces(dart));
        currentLevel_.cellStatistics.postMergeFaces(dart, result);
        return result;
    }

    FaceInfo &removeBridgeInternal(const DartTraverser & dart)
    {
        currentLevel_.cellStatistics.preRemoveBridge(dart);
        FaceInfo &result(currentLevel_.segmentation.removeBridge(dart));
        currentLevel_.cellStatistics.postRemoveBridge(dart, result);
        return result;
    }

    EdgeInfo &mergeEdgesInternal(const DartTraverser & dart)
    {
        currentLevel_.cellStatistics.preMergeEdges(dart);
        EdgeInfo &result(currentLevel_.segmentation.mergeEdges(dart));
        currentLevel_.cellStatistics.postMergeEdges(dart, result);
        return result;
    }

    CellInfo &performOperation(Operation &op,
                               Level *levelData)
    {
        DartTraverser param(&levelData->segmentation, op.param);

        switch(op.type)
        {
          case RemoveIsolatedNode:
          {
              return removeIsolatedNodeInternal(param);
          }
          case MergeFaces:
          {
              return mergeFacesInternal(param);
          }
          case RemoveBridge:
          {
              return removeBridgeInternal(param);
          }
          case MergeEdges:
          {
              return mergeEdgesInternal(param);
          }
          case RemoveEdge:
          {
              return (param.leftFaceLabel() == param.rightFaceLabel() ?
                      removeBridgeInternal(param) :
                      mergeFacesInternal(param));
          }
          case RemoveEdgeWithEnds:
          {
              EdgeInfo &removedEdge = levelData->segmentation.edge(param.edgeLabel());
              NodeInfo &node1(removedEdge.start.startNode());
              NodeInfo &node2(removedEdge.end.startNode());

              FaceInfo &result = (param.leftFaceLabel() == param.rightFaceLabel() ?
                                  removeBridgeInternal(param) :
                                  mergeFacesInternal(param));

              if(!node1.degree)
                  removeIsolatedNodeInternal(node1.anchor);
              /*else
              {
                  DartTraverser changedNode(node1.anchor);
                  // test if the node has degree 2 and the edge is no
                  // loop:
                  if((changedNode.nextSigma() != node1.anchor) &&
                     (changedNode.edgeLabel() != node1.anchor.edgeLabel()) &&
                     (changedNode.nextSigma() == node1.anchor))
                  {
                      mergeEdgesInternal(changedNode);
                  }
                  }*/

              bool removedEdgeIsLoop = (node1.label == node2.label);
              if(!removedEdgeIsLoop)
              {
                  if(!node2.degree)
                      removeIsolatedNodeInternal(node2.anchor);
                  /*else
                  {
                      DartTraverser changedNode(node2.anchor);
                      if((changedNode.nextSigma() != node2.anchor) &&
                         (changedNode.edgeLabel() != node2.anchor.edgeLabel()) &&
                         (changedNode.nextSigma() == node2.anchor))
                      {
                          mergeEdgesInternal(changedNode);
                      }
                      }*/
              }

              return result;
          }
        }

        vigra_fail("Unknown operation type in CellPyramid<>::performOperation!");
        return levelData->segmentation.face(0);
    }

        /** Returns false (and does not change currentLevel()) if
         * currentLevel() is a "better" position for reaching level
         * than the last checkpoint, that is:
         *
         * lastCheckpointIt->first < currentLevel() < level
         */
    bool gotoLastCheckpointBefore(unsigned int level,
                                  Level *levelData)
    {
        vigra_precondition(level < levelCount(),
                           "trying to access non-existing level (too high index)!");

        typename CheckpointMap::iterator lastCheckpointIt =
            checkpoints_.upper_bound(level);
        --lastCheckpointIt;

        if((levelData->level < level) && (lastCheckpointIt->first < levelData->level))
            return false;

        *levelData = lastCheckpointIt->second;
        return true;
    }

    CellInfo &changeLevelInternal(unsigned int gotoLevel,
                                  Level *levelData)
    {
        gotoLastCheckpointBefore(gotoLevel, levelData);

        CellInfo *result = &levelData->segmentation.face(0);

        while(levelData->level < gotoLevel)
        {
            result = &performOperation(history_[levelData->level], levelData);
            ++levelData->level;
        }

        return *result;
    }

  public:
    Level &makeCheckpoint()
    {
        if(!checkpoints_.count(currentLevel()))
            checkpoints_.insert(std::make_pair(currentLevel(), currentLevel_));
        return checkpoints_.find(currentLevel())->second;
    }

    void checkpointChanged(unsigned int which)
    {
        typename CheckpointMap::iterator lastCheckpointIt =
            checkpoints_.upper_bound(currentLevel());
        if((--lastCheckpointIt)->first == which)
        {
            unsigned int oldLevel = currentLevel();
            currentLevel_ = lastCheckpointIt->second;
            changeLevelInternal(oldLevel, &currentLevel_);
        }
    }

    CellPyramid(const Segmentation &level0,
                const CellStatistics &level0Stats = CellStatistics())
    : currentLevel_(0, level0, level0Stats)
    {
        makeCheckpoint();
    }

    FaceInfo &removeIsolatedNode(const DartTraverser & dart)
    {
        history_.push_back(Operation(RemoveIsolatedNode, dart));
        return static_cast<FaceInfo &>(changeLevelInternal(levelCount()-1,
                                                           &currentLevel_));
    }

    FaceInfo &mergeFaces(const DartTraverser & dart)
    {
        history_.push_back(Operation(MergeFaces, dart));
        return static_cast<FaceInfo &>(changeLevelInternal(levelCount()-1,
                                                           &currentLevel_));
    }

    FaceInfo &removeBridge(const DartTraverser & dart)
    {
        history_.push_back(Operation(RemoveBridge, dart));
        return static_cast<FaceInfo &>(changeLevelInternal(levelCount()-1,
                                                           &currentLevel_));
    }

    EdgeInfo &mergeEdges(const DartTraverser & dart)
    {
        history_.push_back(Operation(MergeEdges, dart));
        return static_cast<EdgeInfo &>(changeLevelInternal(levelCount()-1,
                                                           &currentLevel_));
    }

    FaceInfo &removeEdge(const DartTraverser & dart)
    {
        history_.push_back(Operation(RemoveEdge, dart));
        return static_cast<FaceInfo &>(changeLevelInternal(levelCount()-1,
                                                           &currentLevel_));
    }

    FaceInfo &removeEdgeWithEnds(const DartTraverser & dart)
    {
        history_.push_back(Operation(RemoveEdgeWithEnds, dart));
        return static_cast<FaceInfo &>(changeLevelInternal(levelCount()-1,
                                                           &currentLevel_));
    }

    const unsigned int currentLevel() const
    {
        return currentLevel_.level;
    }

    const Segmentation &currentSegmentation() const
    {
        return currentLevel_.segmentation;
    }

    const CellStatistics &currentCellStatistics() const
    {
        return currentLevel_.cellStatistics;
    }

    const Level &gotoLevel(unsigned int level)
    {
        changeLevelInternal(level, &currentLevel_);
        return currentLevel_;
    }

        /** Do a maximum of maxSteps operations to reach given level.
         * Returns true if that was enough, that is (currentLevel() ==
         * level)
         */
    bool approachLevel(unsigned int level, unsigned int maxSteps = 20)
    {
        unsigned int step =
            gotoLastCheckpointBefore(level, &currentLevel_) ? 1 : 0;

        while((step++ < maxSteps) && (currentLevel_.level < level))
        {
            performOperation(history_[currentLevel_.level], &currentLevel_);
            ++currentLevel_.level;
        }

        return (currentLevel_.level == level);
    }

    unsigned int levelCount() const
    {
        return history_.size() + 1;
    }

    void cutHead()
    {
        history_.erase(history_.begin() + currentLevel_.level, history_.end());
        checkpoints_.erase(checkpoints_.upper_bound(currentLevel_.level),
                           checkpoints_.end());
    }
};

} // namespace vigra

#endif // CELLPYRAMID_HXX
