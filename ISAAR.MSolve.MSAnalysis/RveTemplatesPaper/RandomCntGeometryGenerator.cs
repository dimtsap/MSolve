using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using Troschuetz.Random;

namespace ISAAR.MSolve.MSAnalysis.RveTemplatesPaper
{
    /// <summary>
    /// Authors: Odysseas Kokkinos
    /// </summary>
    public class RandomCntGeometryGenerator
    {
        private readonly int _numberOfSimulations;
        private readonly int _numberOfElementsPerCnt;
        private readonly double _cntLength;
        private readonly double _elementLength;
        private readonly int _numberOfCnts;
        private readonly double _standardDeviation;
        private readonly double _upperAngleBound;
        private readonly double _matrixLength;
        private readonly double _matrixWidth;
        private readonly double _matrixHeight;

        public RandomCntGeometryGenerator(int numberOfSimulations,
            int numberOfElementsPerCnt,double cntLength,int numberOfCnts,
            double standardDeviation, double upperAngleBound,
            double matrixLength, double matrixWidth, double matrixHeight)
        {
            _numberOfSimulations = numberOfSimulations;
            _numberOfElementsPerCnt = numberOfElementsPerCnt;
            _cntLength = cntLength;
            if (_numberOfElementsPerCnt <= 0)
                throw new ArgumentOutOfRangeException("The number of CNTS must be greater than zero");

            _elementLength = _cntLength / (double) _numberOfElementsPerCnt;
            _numberOfCnts = numberOfCnts;
            _standardDeviation = standardDeviation;
            _upperAngleBound = Math.Abs(upperAngleBound);
            _matrixLength = matrixLength;
            _matrixWidth = matrixWidth;
            _matrixHeight = matrixHeight;
        }

        public void GenerateCnts()
        {
            var random = new TRandom();
            var numberOfNodesPerCnt = _numberOfElementsPerCnt + 1;
            var nodeCoordinatePerCnt = new double[numberOfNodesPerCnt, 3];
            var cntNodalCoordinates = new double[_numberOfCnts * _numberOfElementsPerCnt, 3];


            for (int indexSimulation = 0; indexSimulation < _numberOfSimulations; indexSimulation++)
            {
                for (int indexCnt = 0; indexCnt < _numberOfCnts; indexCnt++)
                {
                    Array.Clear(nodeCoordinatePerCnt,0,numberOfNodesPerCnt*3);

                    var iNode0 = indexCnt * numberOfNodesPerCnt;
                    cntNodalCoordinates[iNode0,0] = random.ContinuousUniform(0.0, _matrixLength);
                    cntNodalCoordinates[iNode0,1] = random.ContinuousUniform(0.0, _matrixHeight);
                    cntNodalCoordinates[iNode0,2] = random.ContinuousUniform(0.0, _matrixWidth);
                    
                    var x0 = cntNodalCoordinates[iNode0 ,0];
                    var y0 = cntNodalCoordinates[iNode0, 1];
                    var z0 = cntNodalCoordinates[iNode0, 2];

                    for (int indexElement = 0; indexElement < _numberOfElementsPerCnt; indexElement++)
                    {
                        double dx;
                        double dy;
                        double dz;
                        var iNode = indexCnt * numberOfNodesPerCnt + indexElement;
                        if (indexElement == 0)
                        {
                            var thetaXY = random.ContinuousUniform(-Math.PI, Math.PI);
                            var thetaXZ= random.ContinuousUniform(-Math.PI, Math.PI);
                             dx = Math.Cos(thetaXY) * Math.Cos(thetaXZ) * _elementLength;
                             dy = Math.Sin(thetaXY) * Math.Cos(thetaXZ) * _elementLength;
                             dz = Math.Sin(thetaXZ) * _elementLength;
                        }
                        else
                        {
                            var thetaXY = random.Normal(0.0, 1)*_standardDeviation;
                            var thetaXZ = random.Normal(0.0, 1) * _standardDeviation;
                            if (Math.Abs(thetaXY) > _upperAngleBound)
                                thetaXY = Math.Sign(_upperAngleBound);

                            if (Math.Abs(thetaXZ) > _upperAngleBound)
                                thetaXZ = Math.Sign(_upperAngleBound);

                            var phiXY = Math.Atan2(
                                cntNodalCoordinates[indexElement, 1] - cntNodalCoordinates[indexElement - 1, 1],
                                cntNodalCoordinates[indexElement, 0] - cntNodalCoordinates[indexElement - 1, 0]);
                            var phiXZ = Math.Atan2(
                                cntNodalCoordinates[indexElement, 2] - cntNodalCoordinates[indexElement - 1, 2],
                                Math.Sqrt(
                                    Math.Pow(cntNodalCoordinates[indexElement, 0] - cntNodalCoordinates[indexElement - 1, 0], 2) +
                                    Math.Pow(cntNodalCoordinates[indexElement, 1] - cntNodalCoordinates[indexElement - 1, 1], 2)));

                             dx = Math.Cos(thetaXY + phiXY) * Math.Cos(thetaXZ + phiXZ);
                             dy = Math.Sin(thetaXY + phiXY) * Math.Cos(thetaXZ + phiXZ) * _elementLength;
                             dz = Math.Sin(thetaXZ + phiXZ) * _elementLength;
                        }

                        var distance = Math.Sqrt(dx * dx + dy * dy + dz * dz);

                        cntNodalCoordinates[indexElement, 0] += dx;
                        cntNodalCoordinates[indexElement, 1] += dy;
                        cntNodalCoordinates[indexElement, 2] += dz;


                        var xStart = cntNodalCoordinates[indexElement - 1, 0];
                        var yStart = cntNodalCoordinates[indexElement - 1, 1];
                        var zStart = cntNodalCoordinates[indexElement - 1, 2];
                        var xEnd = cntNodalCoordinates[indexElement, 0];
                        var yEnd = cntNodalCoordinates[indexElement, 1];
                        var zEnd = cntNodalCoordinates[indexElement, 2];

                        cntNodalCoordinates[iNode, 0] =
                            cntNodalCoordinates[iNode-1, 0] + cntNodalCoordinates[indexElement, 0];
                        cntNodalCoordinates[iNode, 1] =
                            cntNodalCoordinates[iNode-1, 1] + cntNodalCoordinates[indexElement, 1];
                        cntNodalCoordinates[iNode, 2] =
                            cntNodalCoordinates[iNode-1, 2] + cntNodalCoordinates[indexElement, 2];

                        if (cntNodalCoordinates[iNode, 0] > _matrixLength || cntNodalCoordinates[iNode, 0] < 0.0 
                           || cntNodalCoordinates[iNode, 1] > _matrixHeight || cntNodalCoordinates[iNode, 1] < 0.0
                           || cntNodalCoordinates[iNode, 2] > _matrixWidth || cntNodalCoordinates[iNode, 2] < 0.0)
                        {
                            indexElement = 0; 
                            continue;
                        }
                    }
                }
            }
        }



    }
}
