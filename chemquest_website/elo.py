def rate_single(rating_a: float, rating_b: float, score_a: float, score_b: float, max_adj: float = 32.0) -> tuple[float, float]:
    """
    Performs Elo rating adjustments for a single game.
    Implementation based on https://en.wikipedia.org/wiki/Elo_rating_system#Mathematical_details.

    Args:
        rating_a (float): Rating of player A
        rating_b (float): Rating of player B
        score_a (float): Score of player A (0.0 to 1.0)
        score_b (float): Score of player B (0.0 to 1.0)
        max_adj (float, optional): Maximum rating adjustment. Defaults to 32.0.

    Returns:
        tuple[float, float]: New ratings for player A and player B
    """
    expected_score_a = 1/(1 + 10 ** ((rating_b - rating_a) / 400))
    expected_score_b = 1/(1 + 10 ** ((rating_a - rating_b) / 400))

    new_rating_a = rating_a + max_adj * (score_a - expected_score_a)
    new_rating_b = rating_b + max_adj * (score_b - expected_score_b)

    return (new_rating_a, new_rating_b)

def calculate_new_elo(old_elo: float, question_elo: float, attempts: list[bool]) -> float:
    """
    Calculates the new Elo of a player after attempting a question.

    The player's score is calculated based on the number of attempts and whether the question was ultimately answered correctly.
    A score of 1.0 is awarded for answering correctly on the first attempt, and 0.4 is deducted for each incorrect attempt, with a minimum score of 0.0.
    Attempts after the correct answer are not considered.

    Then, the player's Elo is adjusted based on the question's Elo and the player's score.

    Args:
        old_elo (float): Player's elo before attempting the question
        question_elo (float): Question's elo
        attempts (list[bool]): List of attempts, True for correct, False for incorrect

    Returns:
        float: New elo of the player
    """
    score = 0.0
    for attempt in attempts:
        if attempt == False:
            score -= 0.4
        else:
            score += 1.0
            break
    
    score = max(0.0, score)
    
    new_elo = rate_single(old_elo, question_elo, score, 1.0 - score)[0]
    return new_elo