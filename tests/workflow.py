"""Public, governed-data-free orchestration contract tests."""
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("workflow", REPO / "analysis/workflow/run.py")
w = importlib.util.module_from_spec(spec); spec.loader.exec_module(w)


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.c = json.loads((w.HERE / "manuscript.json").read_text())
        self.c["root"] = self.tmp.name

    def test_graph(self):
        seen = set()
        for stage, parents in w.PARENTS.items():
            self.assertTrue(set(parents) <= seen)
            seen.add(stage)
        self.assertEqual(w.selection(["finalize"]), list(w.PARENTS))
        self.assertIn("paternal", w.selection(["supplementary"]))
        self.assertIn("upper", w.selection(["checklist"]))

    def test_commands_and_manifest_routes(self):
        # Model job paths are produced during prepare; exercise the real reader.
        for mode,n in zip(w.MODES,(195,10,30,10)):
            p = Path(self.tmp.name)/f"jobs/{self.c['run_id']}/{mode}/tasks.txt"
            p.parent.mkdir(parents=True); p.write_text("\n".join(f"/test/{mode}_{i}.rds" for i in range(n)))
        for stage in w.PARENTS:
            for command in w.commands(stage,self.c,1 if stage in w.ARRAYS else None):
                self.assertTrue(Path(command[2]).is_file(),command)
                self.assertNotIn("R2Nov", " ".join(command))
        commands = w.commands("cohort",self.c)[0]
        self.assertIn("--model-manifest",commands)
        self.assertTrue(commands[-1].endswith("model_manifest_main_manuscript_r2.rds"))
        self.assertIn("--cohort-summary-file",w.commands("paternal",self.c)[0])

    def test_identity_and_modified_outputs(self):
        w.initialize(self.c)
        p=Path(self.tmp.name)/"outputs/test.txt"; p.parent.mkdir(); p.write_text("ok")
        record={"identity":w.sha(Path(self.tmp.name)/"workflow/identity.json"),"outputs":{"outputs/test.txt":w.sha(p)}}
        w.atomic_json(w.record_path(self.c,"main"),record)
        w.verify(self.c,"main")
        p.write_text("changed")
        with self.assertRaises(RuntimeError): w.verify(self.c,"main")
        self.c["model_seed"] += 1
        with self.assertRaises(RuntimeError): w.initialize(self.c)

    def test_slurm_dependency_and_array(self):
        w.initialize(self.c); (Path(self.tmp.name)/"logs").mkdir()
        calls=[]
        def sbatch(args,**kwargs):
            calls.append((args,kwargs["input"]))
            return type("Result",(),{"stdout":str(100+len(calls))})()
        with patch.object(w.subprocess,"run",side_effect=sbatch):
            w.submit(self.c,["prepare","model_fit","models"],set())
        self.assertIn("--array=1-245",calls[1][0])
        self.assertIn("--dependency=afterok:101",calls[1][0])
        self.assertIn("--dependency=afterok:102",calls[2][0])
        self.assertIn('--task "$SLURM_ARRAY_TASK_ID"',calls[1][1])
        with self.assertRaises(RuntimeError): w.submit(self.c,["prepare"],set())

    def test_failed_stage_does_not_publish(self):
        with patch.object(w,"commands",return_value=[["false"]]):
            with self.assertRaises(w.subprocess.CalledProcessError): w.execute(self.c,"prepare")
        self.assertFalse(w.record_path(self.c,"prepare").exists())


if __name__ == "__main__": unittest.main()
