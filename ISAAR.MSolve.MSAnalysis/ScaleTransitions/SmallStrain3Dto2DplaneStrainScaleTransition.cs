﻿using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.FEM.Interfaces;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.MultiscaleAnalysis
{
    public class SmallStrain3Dto2DplaneStrainScaleTransition : IScaleTransitions_v2
    {
        public SmallStrain3Dto2DplaneStrainScaleTransition()
        { }

        public double[] MacroToMicroTransition(Node_v2 boundaryNode, double[] MacroScaleVariable)
        {
            double[,] Dq_nodal = new double[3,2]; // Prosoxh: pithanes diorthoseis eis triploun
            Dq_nodal[0, +0] = boundaryNode.X;
            Dq_nodal[1, +1] = boundaryNode.Y;
            Dq_nodal[2, +0] = 0.5*boundaryNode.Y;
            Dq_nodal[2, +1] = 0.5 * boundaryNode.X;


            double[] microVariable = new double[2];            

            for (int i1 = 0; i1 < 2; i1++)
            {
                for (int j1 = 0; j1 < 3; j1++)
                {
                    microVariable[i1] += Dq_nodal[j1, i1] * MacroScaleVariable[j1]; //einai sunolikh 
                }
            }

            return microVariable;
        }

        public double[] MicroToMacroTransition(INode boundaryNode, double[] MicroScaleVariable)
        {
            double[,] Dq_nodal = new double[3, 2];
            Dq_nodal[0, +0] = boundaryNode.X;
            Dq_nodal[1, +1] = boundaryNode.Y;
            Dq_nodal[2, +0] = 0.5 * boundaryNode.Y;
            Dq_nodal[2, +1] = 0.5 * boundaryNode.X;

            double[] macroVariable = new double[3];
            //
            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int j1 = 0; j1 < 2; j1++)
                {
                    macroVariable[i1] += Dq_nodal[ i1, j1] * MicroScaleVariable[j1]; //einai sunolikh 
                }
            }

            return macroVariable;
        }

        public int PrescribedDofsPerNode()
        {
            return 2;
        }

        public int MacroscaleVariableDimension()
        {
            return 3;
        }

        public void ModifyMicrostructureTotalPrescribedBoundaryDisplacementsVectorForMacroStrainVariable(Node_v2 boundaryNode,
            double[] smallStrain2Dmacro, Dictionary<int, Dictionary<DOFType, double>> totalPrescribedBoundaryDisplacements)
        {
            //double[,] Dq_nodal = new double[9, 3];
            double[,] Dq_nodal = new double[3, 2];
            Dq_nodal[0, +0] = boundaryNode.X;
            Dq_nodal[1, +1] = boundaryNode.Y;
            Dq_nodal[2, +0] = 0.5 * boundaryNode.Y;
            Dq_nodal[2, +1] = 0.5 * boundaryNode.X;

            //double[] thesi_prescr_xyz = new double[2];
            double[] u_prescr_xyz_sunol = new double[2];

            for (int i1 = 0; i1 < 2; i1++)
            {
                for (int j1 = 0; j1 < 3; j1++)
                {
                    u_prescr_xyz_sunol[i1] += Dq_nodal[j1, i1] * smallStrain2Dmacro[j1]; //einai sunolikh 
                }
            }

            //SHMEIWSH: an prosthesoume sto totalBoundaryNodalDIsplacements trith metakinhsh (dld logw u_prescr_xyz_sunol[3]) ==0
            // ousiastika h methodos "ImposePrescribedDisplacementsWithInitialConditionSEffect" ths subdomain.cs tha paei kai tha 
            //xanagrapsei panw apo to 0 tou localsolution enos element (pou prokuptei ean exoume thesei sthnepomenh methodo dof = constrained) thn timh 0
            //pou tha vrei sto totalBoundaryNodalDIsplacements. an ekei den exoume valei timh (dld logw u_prescr_xyz_sunol[3] ==0
            //sto  DOFType.Z tou Dictionary) den tha peiraxei to mhden tou localsolution opote einai to idio.

            Dictionary<DOFType, double> totalBoundaryNodalDisplacements = new Dictionary<DOFType, double>();
            totalBoundaryNodalDisplacements.Add(DOFType.X, u_prescr_xyz_sunol[0]);
            totalBoundaryNodalDisplacements.Add(DOFType.Y, u_prescr_xyz_sunol[1]);
            //totalBoundaryNodalDisplacements.Add(DOFType.Z, u_prescr_xyz_sunol[2]); (OPOTE plane strain mporoume na ton valiume ton periorismo apla tha epanaupologizetai)

            totalPrescribedBoundaryDisplacements.Add(boundaryNode.ID, totalBoundaryNodalDisplacements);
        }

        public void ImposeAppropriateAndRigidBodyConstraintsPerBoundaryNode(Model_v2 model, Node_v2 boundaryNode, Dictionary<Node_v2, IList<DOFType>> RigidBodyNodeConstraints)
        {
                model.NodesDictionary[boundaryNode.ID].Constraints.Add(new Constraint() {DOF=DOFType.X});
                model.NodesDictionary[boundaryNode.ID].Constraints.Add(new Constraint() { DOF = DOFType.Y });
                model.NodesDictionary[boundaryNode.ID].Constraints.Add(new Constraint() { DOF = DOFType.Z });           
        }

        public void ImposeAppropriateConstraintsPerBoundaryNode(Model_v2 model, Node_v2 boundaryNode)
        {
            throw new System.NotSupportedException();
        }
    }
}
