import acts
import acts.examples
from AliceActsPythonBindings import ExampleAlgorithm


s = acts.examples.Sequencer(
    events=10, trackFpes=False, numThreads=1, outputDir="./")


cfg = ExampleAlgorithm.Config()
cfg.testString = "your test string value"

customLogLevel = acts.examples.defaultLogging(s, acts.logging.INFO)

exampleAlgorithm = ExampleAlgorithm(config=cfg,
                                    level=acts.logging.INFO)


s.addAlgorithm(exampleAlgorithm)

s.run()
