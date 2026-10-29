library(stemfeedback)

model <- minimal_model(N = 100, p = 1, q = 4, a = .75, dy = 4)
initial <- start_from_stems(25, 25, model)

prediction <- diffusion_probability(initial, model)
print(prediction)

reference <- solve_reference(model, limits = c(150, 150))
result <- simulate_competition(model, initial, n = 100, seed = 123,
                               reference = reference)
print(result)
print(result$feedback)

invasion <- minimal_model(2000, q = 1, b = 1, dy = 1)
print(latency_shift(invasion))

print(recovery_prediction(c(25, 12.5), minimal_model(q = 3, dy = 1.5),
                          exposure = c(3, 3)))
print(neutral_mean_time(.5, minimal_model(q = 1, dy = 1)))
